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
from typing import Dict, List, Tuple, Union, Any
import warnings

# =============================================================================
# External Python modules
# =============================================================================
import numpy as np


# =============================================================================
# Extension modules
# =============================================================================
from . import MExt


class AVLSolver(object):
    def __init__(self, geo_file=None, mass_file=None, debug=False, timing=False):
        # This is important for creating multiple instances of the AVL solver that do not share memory
        curDir = os.path.basename(os.path.dirname(os.path.realpath(__file__)))
        time.sleep(0.1)  # this is necessary for some reason?!
        self.avl = MExt.MExt("libavl", curDir, debug=debug)._module

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
            
            if timing:
                self.avl.case_l.ltiming = True

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

        #  the case parameters are stored in a 1d array,
        # these indices correspond to the position of each parameter in that arra
        self.param_idx_dict = {
            "alpha": 0,
            "beta": 1,
            "pb/2V": 2,
            "qc/2V": 3,
            "rb/2V": 4,
            "CL": 5,
            "CD0": 6,
            "bank": 7,
            "elevation": 8,
            "heading": 9,
            "Mach": 10,
            "velocity": 11,
            "density": 12,
            "grav.acc.": 13,
            "turn rad.": 14,
            "load fac.": 15,
            "X cg": 16,
            "Y cg": 17,
            "Z cg": 18,
            "mass": 19,
            "Ixx": 20,
            "Iyy": 21,
            "Izz": 22,
            "Ixy": 23,
            "Iyz": 24,
            "Izx": 25,
            "visc CL_a": 26,
            "visc CL_u": 27,
            "visc CM_a": 28,
            "visc CM_u": 29,
        }

    def add_constraint(self, var, val, con_var=None):
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
        options = {
            "bankAng": ["B"],
            "CL": ["C"],
            "velocity": ["V"],
            "mass": ["M"],
            "dens": ["D"],
            "G": ["G"],
            "X cg": ["X"],
            "Y cg": ["Y"],
            "Z cg": ["Z"],
        }

        if not (variable in options):
            raise ValueError(
                f"constraint variable `{variable}` not a valid option. Must be one of the following {[key for key in options]} "
            )

        self.avl.trmset("C1", "1 ", options[variable][0], (str(val) + "  \n"))

    def get_case_total_data(self) -> Dict[str, float]:
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
            "CM": ["CASE_R", "CMTOT"],
            "CN BA": ["CASE_R", "CNTOT"],

            # non-dimensionalized moments (stablity frame)
            "CR SA": ["CASE_R", "CRSAX"],
            # "CM SA": ["CASE_R", "CMSAX"], # This is the same in both frames
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

    def set_avl_fort_var(self, common_block, variable, val):
        # convert from fortran ordering to c ordering
        val = val.ravel(order="C").reshape(val.shape[::-1], order="F")

        # this had to be split up into two steps to work
        # get the corresponding common block object.
        # it must be lowercase becuase of f2py
        common_block = getattr(self.avl, common_block.lower())

        # get the value of the variable from the common block
        setattr(common_block, variable.lower(), val)

        return val

    def get_case_surface_data(self) -> Dict[str, Dict[str, float]]:
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

    def get_case_parameter(self, param_key: str) -> float:
        """
        analogous to ruinont Modify parameters for the oper menu to view parameters.
        """
        parvals = self.get_avl_fort_var("CASE_R", "PARVAL")

        # [0] becuase pyavl only supports 1 run case
        param_val = parvals[0][self.param_idx_dict[param_key]]

        return param_val

    def set_case_parameter(self, param_key: str, param_val: float) -> None:
        """
        analogous to ruinont Modify parameters for the oper menu to view parameters.
        """
        # warn the user that alpha, beta,
        if param_key is ["alpha", "beta", "pb/2V", "qc/2V", "rb/2V", "CL"]:
            raise ValueError(
                "alpha, beta, pb/2V, qc/2V, rb/2V, and CL are not allowed to be set,\n\
                             they are calculated during each run based on the contraints. to specify\n\
                             one of these values use the set_contrain method."
            )

        parvals = self.get_avl_fort_var("CASE_R", "PARVAL")
        # [0] becuase pyavl only supports 1 run case
        parvals[0][self.param_idx_dict[param_key]] = param_val

        self.set_avl_fort_var("CASE_R", "PARVAL", parvals)

        # (1) here becuase we want to set the first runcase with fortran indexing (the only one)
        self.avl.set_params(1)

    # def get_condition_data(self) -> Dict[str, float]:
    #     """get the data defining the run condition from the fortran layer"""
    #     # fmt: off
    #     var_to_avl_var = {
    #         "alpha": ["CASE_R", "ALPHA"],
    #         "beta": ["CASE_R", "BETA"],
    #         "Vinf": ["CASE_R", "VINF"],
    #         "mach": ["CASE_R", "MACH"],
    #     }
    #     # fmt: on
    #     condition_data = {}

    #     for key, avl_key in var_to_avl_var.items():
    #         val = self.get_avl_fort_var(*avl_key)
    #         # [()] becuase all the data is stored as a ndarray.
    #         # for scalars this results in a 0-d array.
    #         # It is easier to work with floats so we extract the value with [()]
    #         condition_data[key] = val[()]

    #     return condition_data

    def get_strip_data(self) -> Dict[str, Dict[str, np.ndarray]]:
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
        warnings.warn("executeRun is deprecated, use execute_run instead")
        self.execute_run()

    def execute_run(self):
        # run the analysis (equivalent to the avl command `x` in the oper menu)
        self.avl.oper()

        # needed for stability derivaties
        self.avl.calcst()

    def CLSweep(self, start_CL, end_CL, increment=0.1):
        CLs = np.arange(start_CL, end_CL + increment, increment)

        for cl in CLs:
            self.add_trim_condition("CL", cl)
            self.execute_run()

    def get_control_names(self) -> List[str]:
        control_names = self._convertFortranStringArrayToList(self.avl.case_c.dname)
        return control_names

    def get_surface_names(self) -> List[str]:
        """get the surface names from the geometry"""
        surf_names = self._convertFortranStringArrayToList(self.avl.case_c.stitle)
        return surf_names

    def get_surface_params(self, geom_only: bool = False) -> Dict[str, Dict[str, Any]]:
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

    def set_surface_params(self, surf_data: Dict[str, Dict[str, any]]) -> None:
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

    # Utility functions
    def get_num_surfaces(self) -> int:
        """Get the number of surfaces in the geometry"""
        return self.avl.case_i.nsurf

    def get_num_strips(self) -> int:
        """Get the number of strips in the mesh"""
        return self.avl.case_i.nstrip

    def get_mesh_size(self) -> int:
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

    def solve():
        """
            use scipy's linear solves to solve the system of equations. 
        """
        pass
        

class CaseData:
    """class to hold the resulting data from a single run.

    The data is stored in dictionaries that are indexed by the surface name.
    The dictionary names are related to oper command used to show the data.

    """

    def __init__(self) -> None:
        # Data used to normalize the other data (Sref, Cref, Bref, etc.)
        self.reference_data = {}

        # Data of the forces and moments for everything (total_data), each surface (surface_data), each body (body_data), each strip (strip_data), and each element (element_data)
        self.total_data = {}
        self.surface_data = {}
        self.body_data = {}
        self.strip_data = {}
        self.element_data = {}

        # Data of the stability derivatives in the stability axis (self.stability_axis_derivs) and body axis (self.body_axis_derivs)
        self.stability_axis_derivs = {}
        self.body_axis_derivs = {}
