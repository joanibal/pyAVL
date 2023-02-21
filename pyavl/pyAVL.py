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

        # store the control surfaces variables
        control_names = self._convertFortranStringArrayToList(self.avl.case_c.dname)

        self.control_variables = {}
        for idx_c_var, c_name in enumerate(control_names):
            self.control_variables[c_name] = f"D{idx_c_var+1}"

    def addConstraint(self, var, val, con_var=None):
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

        self.avl.conset(avl_var, f"{avl_con_var} {val} \n")

    def addTrimCondition(self, variable, val):
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

    def executeRun(self):
        self.__exe = True

        self.avl.oper()

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

        if self.surf_CL.size == 0:
            self.surf_CL = surf_CLs
            self.surf_CD = surf_CDs
        else:
            self.surf_CL = surf_CLs = np.hstack((self.surf_CL, surf_CLs))
            self.surf_CD = surf_CDs = np.hstack((self.surf_CD, surf_CDs))

    def calcNP(self):
        # executeRun must be run first

        if not (self.__exe):
            raise RuntimeError(":  executeRun() must be called first")

        self.avl.calcst()
        self.NP = self.avl.case_r.xnp
        # print 'Xnp:', self.avl.case_r.xnp

    def alphaSweep(self, start_alpha, end_alpha, increment=1):
        alphas = np.arange(start_alpha, end_alpha + increment, increment)

        for alf in alphas:
            self.addConstraint("alpha", alf)
            self.executeRun()

    def CLSweep(self, start_CL, end_CL, increment=0.1):
        CLs = np.arange(start_CL, end_CL + increment, increment)

        for cl in CLs:
            self.addTrimCondition("CL", cl)
            self.executeRun()

    def get_surface_names(self) -> list[str]:
        """get the surface names from the geometry"""
        surf_names = self._convertFortranStringArrayToList(self.avl.case_c.stitle)
        return surf_names

    def get_surface_params(self) -> dict[str, dict[str:any]]:
        """get the surface level parameters from the geometry"""
        surf_names = self.get_surface_names()
        surf_data = {}
        for idx_surf, surf_name in enumerate(surf_names):
            # only copy unduplicated sufaces
            # print(surf_name, self.avl.surf_geom_i.nsec[idx_surf])
            if (self.avl.surf_i.imags[idx_surf] < 0):
                # this is a duplicated surface, skip it
                continue
            
            
            num_sec =  self.avl.surf_geom_i.nsec[idx_surf]

            surf_data[surf_name] = {
                "scale": self.avl.surf_geom_r.xyzscal[:, idx_surf].T,
                "translate": self.avl.surf_geom_r.xyztran[:, idx_surf].T,
                "angle": np.rad2deg(self.avl.surf_geom_r.addinc[idx_surf]),
                "nchordwise":  self.avl.surf_geom_i.nvc[idx_surf],
                "cspace": self.avl.surf_geom_r.cspace[idx_surf],
                "nspan": self.avl.surf_geom_i.nvs[idx_surf],
                "sspace": self.avl.surf_geom_r.sspace[idx_surf],
            }
            if self.avl.surf_geom_l.ldupl[idx_surf]:
                surf_data[surf_name]["yduplicate"] = self.avl.surf_geom_r.ydupl[idx_surf]

        return surf_data

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

    # def calcSectionalCPDist(self):

    #     for N in [1, self.avl.case_i.nsurf]:

    #         J1 = self.avl.surf_i.jfrst(N)
    #         JN = J1 + self.avl.surf_i.nj - 1
    #         JINC = 1

    #         CPSCL = self.avl.   cpfac*  self.avl.case_r.cref

    #         for

    #


#### fortran code

#       DO N = 1, NSURF

#         J1 = JFRST(N)
#         JN = JFRST(N) + NJ(N)-1
#         JINC = 1


#            CPSCL = CPFAC*CREF


#          IP = 0
#          DO J = J1, JN, JINC
#            I1 = IJFRST(J)
#            NV = NVSTRP(J)
#            DO II = 1, NV
#              IV = I1 + II-1
#              XAVE = RV(1,IV)
#              YAVE = RV(2,IV)
#              ZAVE = RV(3,IV)
#              DELYZ = DCP(IV) * CPSCL
#              XLOAD = XAVE
#              YLOAD = YAVE + DELYZ*ENSY(J)
#              ZLOAD = ZAVE + DELYZ*ENSZ(J)

# c             XLOAD = XAVE + DELYZ*ENC(1,IV)
# c             YLOAD = YAVE + DELYZ*ENC(2,IV)
# c             ZLOAD = ZAVE + DELYZ*ENC(3,IV)

#              IF(II.GT.1) THEN
#                IP = IP+1
#                PTS_LINES(1,1,IP) = XLOADOLD
#                PTS_LINES(2,1,IP) = YLOADOLD
#                PTS_LINES(3,1,IP) = ZLOADOLD
#                PTS_LINES(1,2,IP) = XLOAD
#                PTS_LINES(2,2,IP) = YLOAD
#                PTS_LINES(3,2,IP) = ZLOAD
#                ID_LINES(IP) = 0
#              ENDIF
#              XLOADOLD = XLOAD
#              YLOADOLD = YLOAD
#              ZLOADOLD = ZLOAD
#            END DO
#          END DO
#          NLINES = IP
#          NPROJ = 2*NLINES
#          CALL VIEWPROJ(PTS_LINES,NPROJ,PTS_LPROJ)
#          CALL PLOTLINES(NLINES,PTS_LPROJ,ID_LINES)
#          CALL NEWCOLOR(ICOL)
#          CALL NEWPEN(IPN)
#         ENDIF
# C
