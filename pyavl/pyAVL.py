"""
pyAVL

pyAVL is a wrapper for Mark Drela's AVL code. The purpose of this
class is to provide an easy to use wrapper for avl for intergration
into other projects.

Developers:
-----------
- Josh Anibal (JLA)

History
-------
    v. 0.0 - Initial Class Creation (JLA, 08 2016)
    v. 0.1 - Added CDFF and sectional properties for all surfaces (JLA, 12 2017)
    v. 1.0 - major refactor (I learned a lot more programming in 7 years) (JLA, 02 2023)
"""

# =============================================================================
# Standard Python modules
# =============================================================================
import os
import time
import copy
from pprint import pprint

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np


# =============================================================================
# Extension modules
# =============================================================================
from . import MExt


class AVLSolver(object):
    def __init__(self, geo_file=None, mass_file=None, debug=False):
        # This is important for creating multiple instances of the AVL solver that do not share memory
        curDir = os.path.basename(os.path.dirname(os.path.realpath(__file__)))
        time.sleep(0.1)  # this is necessary for some reason?!
        self.avl = MExt.MExt("libavl", curDir, debug=debug)._module

        self.resetData()

        if not (geo_file is None):
            try:
                # check to make sure files exist
                file = geo_file
                f = open(geo_file, "r")
                f.close()

                if not (mass_file is None):
                    file = mass_file
                    f = open(mass_file, "r")
                    f.close()
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"Could not open the file {file} from python. This is usually an issue with the specified file path"
                )

            self.avl.avl()
            if debug:
                self.avl.case_l.lverbose = True

            self.avl.loadgeo(geo_file)

            if mass_file is not None:
                self.avl.loadmass(mass_file)

        else:
            raise ValueError("neither a geometry file or aircraft object was given")

        self.bref = self.avl.case_r.bref
        self.cref = self.avl.case_r.cref
        self.sref = self.avl.case_r.sref

        control_names = self.get_control_names()
        self.control_variables = {}
        for idx_c_var, c_name in enumerate(control_names):
            self.control_variables[c_name] = f"D{idx_c_var+1}"

    def add_constraint(self, var, val, con_var=None):
        self.__exe = False

        avl_variables = {
            "alpha": "A",
            "beta": "B",
            "roll rate": "R",
            "pitch rate": "P",
            "yaw rate": "Y",
        }

        avl_con_variables = copy.deepcopy(avl_variables)
        avl_con_variables.update(
            {
                "CL": "C",
                "CY": "S",
                "Cl roll moment": "RM",
                "Cm pitch moment": "PM",
                "Cn yaw moment": "YM",
            }
        )

        if var in avl_variables:
            # save the name of the avl_var
            avl_var = avl_variables[var]
        elif var in self.control_variables.keys():
            avl_var = self.control_variables[var]
        elif var in self.control_variables.values():
            avl_var = var
        else:
            raise ValueError(
                f"specified variable `{var}` not a valid option. Must be one of the following variables{[key for key in avl_variables]} or control surface name or index{[item for item in self.control_variables.items()]}."
            )

        if con_var is None:
            avl_con_var = avl_var
        elif con_var in avl_con_variables:
            avl_con_var = avl_con_variables[con_var]
        elif con_var in self.control_variables.keys():
            avl_con_var = self.control_variables[con_var]
        elif con_var in self.control_variables.values():
            avl_con_var = con_var
        else:
            raise ValueError(
                f"specified contraint variable `{con_var}` not a valid option. Must be one of the following variables{[key for key in avl_variables]} or control surface name or index{[item for item in self.control_variables.items()]}."
            )

        # check that the type of val is correct
        if not isinstance(val, (int, float, np.floating, np.integer)):
            raise TypeError(f"contraint `val` must be a int or float. Got {type(val)}")

        self.avl.conset(avl_var, f"{avl_con_var} {val} \n")

    def add_trim_condition(self, variable, val):
        self.__exe = False

        options = {
            "bankAng": ["B"],
            "CL": ["C"],
            "velocity": ["V"],
            "mass": ["M"],
            "dens": ["D"],
            "G": ["G"],
            "X_cg": ["X"],
            "Y_cg": ["Y"],
            "Z_cg": ["Z"],
        }

        if not (variable in options):
            raise ValueError(
                f"constraint variable `{variable}` not a valid option. Must be one of the following {[key for key in options]} "
            )

        self.avl.trmset("C1", "1 ", options[variable][0], (str(val) + "  \n"))

    def get_case_total_data(self):
        # fmt: off
        # This dict has the following structure:
        # python key: [common block name, fortran varaiable name]
        var_to_avl_var = {
            # lift and drag from surface integration (wind frame)
            "CL": ["CASE_R", "CLTOT"],
            "CD": ["CASE_R", "CDTOT"],
            "CDv": ["CASE_R", "CDVTOT"], # viscous drag
            
            # lift and drag calculated from farfield integration
            "CLff": ["CASE_R", "CLFF"],
            "CYff": ["CASE_R", "CDFF"],
            "CDi": ["CASE_R", "CDFF"], # viscous drag
            
            # non-dimensionalized forces 
            "CX": ["CASE_R", "CXTOT"], 
            "CY": ["CASE_R", "CYTOT"], 
            "CZ": ["CASE_R", "CZTOT"],   
            
            # non-dimensionalized moments (body frame)
            "CR BA": ["CASE_R", "CRTOT"],
            "CM BA": ["CASE_R", "CMTOT"],
            "CN BA": ["CASE_R", "CNTOT"],

            # non-dimensionalized moments (stablity frame)
            "CR SA": ["CASE_R", "CRSAX"],
            "CM SA": ["CASE_R", "CMSAX"],
            "CN SA": ["CASE_R", "CNSAX"],
            
            # spanwise efficiency
            "e": ["CASE_R", "SPANEF"],
        }
        # fmt: on

        total_data = {}

        for key, avl_key in var_to_avl_var.items():
            val = self.get_avl_fort_var(*avl_key)
            # [()] becuase all the data is stored as a ndarray.
            # for scalars this results in a 0-d array.
            # It is easier to work with floats so we extract the value with [()]
            total_data[key] = val[()]

        return total_data

    def get_avl_fort_var(self, common_block, variable):
        # this had to be split up into two steps to work

        # get the corresponding common block object.
        # it must be lowercase becuase of f2py
        common_block = getattr(self.avl, common_block.lower())

        # get the value of the variable from the common block
        val = getattr(common_block, variable.lower())

        # convert from fortran ordering to c ordering
        val = val.ravel(order="F").reshape(val.shape[::-1], order="C")

        return val

    def get_case_surface_data(self):
        # fmt: off
        # This dict has the following structure:
        # python key: [common block name, fortran varaiable name]
        var_to_avl_var = {
            # surface contributions to total lift and drag from surface integration (wind frame)
            "CL": ["SURF_R", "CLSURF"],
            "CD": ["SURF_R", "CDSURF"],
            "CDv": ["SURF_R", "CDVSURF"],  # viscous drag
            
            # non-dimensionalized forces 
            "CX": ["SURF_R", "CXSURF"], 
            "CY": ["SURF_R", "CYSURF"], 
            "CZ": ["SURF_R", "CZSURF"],   
            
            # non-dimensionalized moments (body frame)
            "CR": ["SURF_R", "CRSURF"],
            "CM": ["SURF_R", "CMSURF"],
            "CN": ["SURF_R", "CNSURF"],
            
            # forces non-dimentionalized by surface quantities
            # uses surface area instead of sref and takes moments about leading edge
            "CL surf" : ["SURF_R", "CL_SRF"],
            "CD surf" : ["SURF_R", "CD_SRF"],
            "CMLE surf" : ["SURF_R", "CMLE_SRF"],
            
            #TODO: add CF_SRF(3,NFMAX), CM_SRF(3,NFMAX)
        }
        # fmt: on

        surf_names = self.get_surface_names()

        # add a dictionary for each surface that will be filled later
        surf_data = {}
        for surf in surf_names:
            surf_data[surf] = {}

        for key, avl_key in var_to_avl_var.items():
            vals = self.get_avl_fort_var(*avl_key)

            # add the values to corresponding surface dict
            for idx_surf, surf_name in enumerate(surf_names):
                surf_data[surf_name][key] = vals[idx_surf]

        return surf_data

    def get_strip_data(self):
        # fmt: off
        var_to_avl_var = {
            # geometric quantities
            "chord": ["STRP_R", "CHORD"],
            "width": ["STRP_R", "WSTRIP"],
            "XYZ LE": ["STRP_R", "RLE"], #control point leading edge coordinates
            "twist": ["STRP_R", "AINC"],
            
            # strip contributions to total lift and drag from strip integration
            "CL": ["STRP_R", "CLSTRP"],
            "CD": ["STRP_R", "CDSTRP"],
            
            # strip contributions to non-dimensionalized forces 
            "CX": ["STRP_R", "CXSTRP"], 
            "CY": ["STRP_R", "CYSTRP"], 
            "CZ": ["STRP_R", "CZSTRP"],   
            
            
            # strip contributions to total moments (body frame)
            "CM": ["STRP_R", "CMSTRP"],
            "CN": ["STRP_R", "CNSTRP"],
            "CR": ["STRP_R", "CRSTRP"],
            
            
            # forces non-dimentionalized by strip quantities
            "CL strip" : ["STRP_R", "CL_LSTRP"],
            "CD strip" : ["STRP_R", "CD_LSTRP"],
            "CF strip" : ["STRP_R", "CF_STRP"], # forces in 3 directions
            "CM strip" : ["STRP_R", "CF_STRP"], # moments in 3 directions
        }    
        
        surf_names = self.get_surface_names()

        # add a dictionary for each surface that will be filled later
        strip_data = {}
        for surf in surf_names:
            strip_data[surf] = {}
        
        
        for key, avl_key in var_to_avl_var.items():
            
            vals = self.get_avl_fort_var(*avl_key)
            
            # add the values to corresponding surface dict
            for idx_surf, surf_name in enumerate(surf_names):
                    
                idx_srp_beg, idx_srp_end = self.get_surface_strip_indices(idx_surf)
                strip_data[surf_name][key] = vals[idx_srp_beg:idx_srp_end]
            
        return strip_data

    def get_surface_strip_indices(self, idx_surf):
        num_strips = np.trim_zeros(self.avl.surf_i.nj)
        idx_srp_beg = np.sum(num_strips[:idx_surf])
        idx_srp_end = np.sum(num_strips[: idx_surf + 1])

        return idx_srp_beg, idx_srp_end

    def executeRun(self):
        self.__exe = True

        # run the analysis (equivalent to the avl command `x` in the oper menu)
        self.avl.oper()
        self.avl.calcst()
        self.get_case_total_data()
        # self.get_case_surface_data()
        # self.get_strip_data()
        # extract the resulting case data form the fortran layer
        # The extracted data is stored in the CaseData object

        self.alpha = np.append(self.alpha, float(self.avl.case_r.alfa))  # *(180.0/np.pi) # returend in radians)
        self.beta = np.append(self.beta, float(self.avl.case_r.beta))  # *(180.0/np.pi) # returend in radians)
        self.CL = np.append(self.CL, float(self.avl.case_r.cltot))
        self.CD = np.append(
            self.CD, float(self.avl.case_r.cdtot)
        )  # = np append(self.avl.case_r.cdvtot)  for total viscous)
        self.CDV = np.append(
            self.CDV, float(self.avl.case_r.cdvtot)
        )  # = np append(self.avl.case_r.cdvtot)  for total viscous)
        self.CDFF = np.append(self.CDFF, float(self.avl.case_r.cdff))
        self.CM = np.append(self.CM, float(self.avl.case_r.cmtot))
        self.span_eff = np.append(self.span_eff, float(self.avl.case_r.spanef))

        deflecs = (np.trim_zeros(self.avl.case_r.delcon)).reshape(len(np.trim_zeros(self.avl.case_r.delcon)), 1)
        # print deflecs,  self.avl.case_r.delcon
        # x.astype(int)

        if self.control_deflection.size == 0:
            self.control_deflection = deflecs
        else:
            self.control_deflection = np.hstack((self.control_deflection, deflecs))

        self.velocity = np.append(self.velocity, np.asarray(self.avl.case_r.vinf))

        # get section properties
        # NS = self.avl.surf_i.nj[0]
        num_strips = np.trim_zeros(self.avl.surf_i.nj)

        # print NS
        # print np.trim_zeros(self.avl.strp_r.clstrp[:])

        # start = 0
        sec_CLs = []
        sec_CDs = []
        sec_CMs = []
        sec_CNs = []
        sec_CRs = []
        sec_Chords = []
        sec_Yles = []
        sec_widths = []
        end = 0
        for i in range(0, len(num_strips)):
            start = end
            end = start + num_strips[i]
            sec_CLs.append(self.avl.strp_r.clstrp[start:end])
            sec_CDs.append(self.avl.strp_r.cdstrp[start:end])
            sec_CMs.append(self.avl.strp_r.cmstrp[start:end])
            sec_CNs.append(self.avl.strp_r.cnstrp[start:end])
            sec_CRs.append(self.avl.strp_r.crstrp[start:end])
            sec_Chords.append(self.avl.strp_r.chord[start:end])
            sec_Yles.append(self.avl.strp_r.rle[1][start:end])
            sec_widths.append(self.avl.strp_r.wstrip[start:end])

        self.sec_CL.append(sec_CLs)
        self.sec_CD.append(sec_CDs)
        self.sec_CM.append(sec_CMs)
        self.sec_CN.append(sec_CNs)
        self.sec_CR.append(sec_CRs)
        self.sec_Chord = sec_Chords
        self.sec_Yle = sec_Yles
        self.sec_width = sec_widths
        surf_CLs = (np.trim_zeros(self.avl.surf_r.clsurf)).reshape(len(np.trim_zeros(self.avl.surf_r.clsurf)), 1)
        surf_CDs = (np.trim_zeros(self.avl.surf_r.cdsurf)).reshape(len(np.trim_zeros(self.avl.surf_r.cdsurf)), 1)

        # if self.surf_CL.size == 0:
        #     self.surf_CL = surf_CLs
        #     self.surf_CD = surf_CDs
        # else:
        #     self.surf_CL = surf_CLs = np.hstack((self.surf_CL, surf_CLs))
        #     self.surf_CD = surf_CDs = np.hstack((self.surf_CD, surf_CDs))

    def calcNP(self):
        # executeRun must be run first

        if not (self.__exe):
            raise RuntimeError(":  executeRun() must be called first")

        self.avl.calcst()
        self.NP = self.avl.case_r.xnp
        # print 'Xnp:', self.avl.case_r.xnp

    # def alphaSweep(self, start_alpha, end_alpha, increment=1):
    #     alphas = np.arange(start_alpha, end_alpha + increment, increment)

    #     case_data = []
    #     for alf in alphas:
    #         self.add_constraint("alpha", alf)
    #         self.executeRun()
    #         case_data.append(self.get_case_total_data())

    def CLSweep(self, start_CL, end_CL, increment=0.1):
        CLs = np.arange(start_CL, end_CL + increment, increment)

        for cl in CLs:
            self.add_trim_condition("CL", cl)
            self.executeRun()

    def get_control_names(self):
        control_names = self._convertFortranStringArrayToList(self.avl.case_c.dname)
        return control_names

    def get_surface_names(self):
        """get the surface names from the geometry"""
        surf_names = self._convertFortranStringArrayToList(self.avl.case_c.stitle)
        return surf_names

    def get_surface_params(self, geom_only: bool = False):
        """get the surface level parameters from the geometry

        geom_only: only return the geometry parameters of the surface
            - xle, yle, zle, chord, twist


        """

        surf_names = self.get_surface_names()
        surf_data = {}
        for idx_surf, surf_name in enumerate(surf_names):
            # only copy unduplicated sufaces
            # print(surf_name, self.avl.surf_geom_i.nsec[idx_surf])
            if self.avl.surf_i.imags[idx_surf] < 0:
                # this is a duplicated surface, skip it
                continue

            num_sec = self.avl.surf_geom_i.nsec[idx_surf]

            # I don't see a case in the code where nasec is not the same for all sections
            # but just in case, we'll do a test and throw an error if not
            nasec = self.avl.surf_geom_i.nasec[:num_sec, idx_surf]
            if np.unique(nasec).size != 1:
                raise RuntimeError("nasec is not the same for all sections")

            nasec = nasec[0]

            icontd = []
            idestd = []
            xhinged = []
            vhinged = []
            gaind = []
            refld = []
            gaing = []
            for idx_sec in range(num_sec):
                nscon = self.avl.surf_geom_i.nscon[idx_sec, idx_surf]
                icontd.append(self.avl.surf_geom_i.icontd[:nscon, idx_sec, idx_surf])
                xhinged.append(self.avl.surf_geom_r.xhinged[:nscon, idx_sec, idx_surf])
                vhinged.append(self.avl.surf_geom_r.vhinged[:, :nscon, idx_sec, idx_surf].T)
                gaind.append(self.avl.surf_geom_r.gaind[:nscon, idx_sec, idx_surf])
                refld.append(self.avl.surf_geom_r.refld[:nscon, idx_sec, idx_surf])

                nsdes = self.avl.surf_geom_i.nsdes[idx_sec, idx_surf]
                idestd.append(self.avl.surf_geom_i.idestd[:nsdes, idx_sec, idx_surf])
                gaing.append(self.avl.surf_geom_r.gaing[:nsdes, idx_sec, idx_surf])

            surf_data[surf_name] = {
                "scale": self.avl.surf_geom_r.xyzscal[:, idx_surf].T,
                "translate": self.avl.surf_geom_r.xyztran[:, idx_surf].T,
                "angle": np.rad2deg(self.avl.surf_geom_r.addinc[idx_surf]),
                "xyzles": self.avl.surf_geom_r.xyzles[:, :num_sec, idx_surf].T,
                "chords": self.avl.surf_geom_r.chords[:num_sec, idx_surf],
                "aincs": self.avl.surf_geom_r.aincs[:num_sec, idx_surf],
                "xasec": self.avl.surf_geom_r.xasec[:nasec, :num_sec, idx_surf].T,
                "sasec": self.avl.surf_geom_r.sasec[:nasec, :num_sec, idx_surf].T,
                "tasec": self.avl.surf_geom_r.tasec[:nasec, :num_sec, idx_surf].T,
            }

            if not geom_only:
                surf_data[surf_name].update(
                    {
                        # paneling parameters
                        "nchordwise": self.avl.surf_geom_i.nvc[idx_surf],
                        "cspace": self.avl.surf_geom_r.cspace[idx_surf],
                        "nspan": self.avl.surf_geom_i.nvs[idx_surf],
                        "sspace": self.avl.surf_geom_r.sspace[idx_surf],
                        "sspaces": self.avl.surf_geom_r.sspaces[:num_sec, idx_surf],
                        "nspans": self.avl.surf_geom_i.nspans[:num_sec, idx_surf],
                        # control surface variables
                        "icontd": icontd,
                        "idestd": idestd,
                        "clcdsec": self.avl.surf_geom_r.clcdsec[:, :num_sec, idx_surf].T,
                        "claf": self.avl.surf_geom_r.claf[:num_sec, idx_surf],
                        "xhinged": xhinged,
                        "vhinged": vhinged,
                        "gaind": gaind,
                        "refld": refld,
                        "gaing": gaing,
                    }
                )

                if self.avl.surf_geom_l.ldupl[idx_surf]:
                    surf_data[surf_name]["yduplicate"] = self.avl.surf_geom_r.ydupl[idx_surf]

        return surf_data

    def set_surface_params(self, surf_data: dict[str, dict[str:any]]):
        """set the give surface data into the current avl object.
        ASSUMES THE CONTROL SURFACE DATA STAYS AT THE SAME LOCATION"""
        surf_names = self.get_surface_names()

        for surf_name in surf_data:
            idx_surf = surf_names.index(surf_name)

            # only set unduplicated sufaces
            if self.avl.surf_i.imags[idx_surf] < 0:
                # this is a duplicated surface, skip it
                continue

            num_sec = self.avl.surf_geom_i.nsec[idx_surf]

            # I don't see a case in the code where nasec is not the same for all sections
            # but just in case, we'll do a test and throw an error if not
            nasec = self.avl.surf_geom_i.nasec[:num_sec, idx_surf]
            if np.unique(nasec).size != 1:
                raise RuntimeError("nasec is not the same for all sections")

            nasec = nasec[0]

            icontd = surf_data[surf_name]["icontd"]
            idestd = surf_data[surf_name]["idestd"]
            xhinged = surf_data[surf_name]["xhinged"]
            vhinged = surf_data[surf_name]["vhinged"]
            gaind = surf_data[surf_name]["gaind"]
            refld = surf_data[surf_name]["refld"]
            gaing = surf_data[surf_name]["gaing"]
            for idx_sec in range(num_sec):
                nscon = self.avl.surf_geom_i.nscon[idx_sec, idx_surf]
                self.avl.surf_geom_i.icontd[:nscon, idx_sec, idx_surf] = icontd[idx_sec]
                self.avl.surf_geom_r.xhinged[:nscon, idx_sec, idx_surf] = xhinged[idx_sec]
                self.avl.surf_geom_r.vhinged[:, :nscon, idx_sec, idx_surf] = vhinged[idx_sec].T
                self.avl.surf_geom_r.gaind[:nscon, idx_sec, idx_surf] = gaind[idx_sec]
                self.avl.surf_geom_r.refld[:nscon, idx_sec, idx_surf] = refld[idx_sec]

                nsdes = self.avl.surf_geom_i.nsdes[idx_sec, idx_surf]
                self.avl.surf_geom_i.idestd[:nsdes, idx_sec, idx_surf] = idestd[idx_sec]
                self.avl.surf_geom_r.gaing[:nsdes, idx_sec, idx_surf] = gaing[idx_sec]

            self.avl.surf_geom_r.xyzscal[:, idx_surf] = surf_data[surf_name]["scale"].T
            self.avl.surf_geom_r.xyztran[:, idx_surf] = surf_data[surf_name]["translate"].T
            self.avl.surf_geom_r.addinc[idx_surf] = np.deg2rad(surf_data[surf_name]["angle"])
            self.avl.surf_geom_i.nvc[idx_surf] = surf_data[surf_name]["nchordwise"]
            self.avl.surf_geom_r.cspace[idx_surf] = surf_data[surf_name]["cspace"]
            self.avl.surf_geom_i.nvs[idx_surf] = surf_data[surf_name]["nspan"]
            self.avl.surf_geom_r.sspace[idx_surf] = surf_data[surf_name]["sspace"]
            self.avl.surf_geom_r.sspaces[:num_sec, idx_surf] = surf_data[surf_name]["sspaces"]
            self.avl.surf_geom_i.nspans[:num_sec, idx_surf] = surf_data[surf_name]["nspans"]
            self.avl.surf_geom_r.xyzles[:, :num_sec, idx_surf] = surf_data[surf_name]["xyzles"].T
            self.avl.surf_geom_r.chords[:num_sec, idx_surf] = surf_data[surf_name]["chords"]
            self.avl.surf_geom_r.aincs[:num_sec, idx_surf] = surf_data[surf_name]["aincs"]
            self.avl.surf_geom_r.xasec[:nasec, :num_sec, idx_surf] = surf_data[surf_name]["xasec"].T
            self.avl.surf_geom_r.sasec[:nasec, :num_sec, idx_surf] = surf_data[surf_name]["sasec"].T
            self.avl.surf_geom_r.tasec[:nasec, :num_sec, idx_surf] = surf_data[surf_name]["tasec"].T
            self.avl.surf_geom_r.clcdsec[:, :num_sec, idx_surf] = surf_data[surf_name]["clcdsec"].T
            self.avl.surf_geom_r.claf[:num_sec, idx_surf] = surf_data[surf_name]["claf"]

        self.avl.update_surfaces()

    # def write_geom_file(self, filename):
    #     """write the current avl geometry to a file"""
    #     surf_geom = self.get_surface_params()

    # def __write_surface(self, fid, surf_name, surf_data):
    #     """write a surface to a file"""
    #     fid.write("SURFACE

    def resetData(self):
        self.__exe = False

        self.alpha = np.zeros(0)
        self.beta = np.zeros(0)
        self.CL = np.zeros(0)
        self.CD = np.zeros(0)
        self.CDV = np.zeros(0)
        self.CDFF = np.zeros(0)
        self.CM = np.zeros(0)
        self.span_eff = np.zeros(0)

        self.elev_def = np.zeros(0)
        self.rud_def = np.zeros(0)

        self.velocity = np.zeros(0)

        self.control_deflection = np.empty(0)

        self.sec_CL = []
        self.sec_CD = []
        self.sec_CM = []
        self.sec_CN = []
        self.sec_CR = []
        self.sec_Chord = []
        self.sec_Yle = []

        self.surf_CL = np.empty(0)
        self.surf_CD = np.empty(0)

    # Utility functions
    def get_num_surfaces(self):
        """Get the number of surfaces in the geometry"""
        return self.avl.case_i.nsurf

    def get_num_strips(self):
        """Get the number of strips in the mesh"""
        return self.avl.case_i.nstrip

    def get_mesh_size(self):
        """Get the number of panels in the mesh"""
        return self.avl.case_i.nvor

    def _createFortranStringArray(self, strList, num_max_char):
        """Setting arrays of strings in Fortran can be kinda nasty. This
        takesa list of strings and returns the array"""

        arr = np.zeros((len(strList), num_max_char), dtype="str")
        arr[:] = " "
        for i, s in enumerate(strList):
            for j in range(len(s)):
                arr[i, j] = s[j]

        return arr

    def _convertFortranStringArrayToList(self, fortArray):
        """Undoes the _createFotranStringArray"""
        strList = []
        for ii in range(len(fortArray)):
            if fortArray[ii].dtype == np.dtype("|S0"):
                # there is characters in the sting to add
                continue
            if fortArray[ii].dtype == np.dtype("|S1"):
                strList.append(b"".join(fortArray[ii]).decode().strip())
            elif fortArray[ii].dtype == np.dtype("<U1"):
                strList.append("".join(fortArray[ii]).strip())
            # TODO: need a more general solution for |S<variable> type
            elif fortArray[ii].dtype == np.dtype("|S16"):
                strList.append(fortArray[ii].decode().strip())
            elif fortArray[ii].dtype == np.dtype("|S40"):
                strList.append(fortArray[ii].decode().strip())
            else:
                raise TypeError(f"Unable to convert {fortArray[ii]} of type {fortArray[ii].dtype} to string")

        return strList


class CaseData:
    """class to hold the resulting data from a single run.

    The data is stored in dictionaries that are indexed by the surface name.
    The dictionary names are related to oper command used to show the data.

    """

    def __init__(self):
        # Data used to normalize the other data (Sref, Cref, Bref, etc.)
        self.reference_data = {}

        # Data used to define the run conditions (alpha, beta, etc.)
        self.condition_data = {}

        # Data of the forces and moments for everything (total_data), each surface (surface_data), each body (body_data), each strip (strip_data), and each element (element_data)
        self.total_data = {}
        self.surface_data = {}
        self.body_data = {}
        self.strip_data = {}
        self.element_data = {}

        # Data of the stability derivatives in the stability axis (self.stability_axis_derivs) and body axis (self.body_axis_derivs)
        self.stability_axis_derivs = {}
        self.body_axis_derivs = {}
