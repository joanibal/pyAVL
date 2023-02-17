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
                    f.close
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"Could not open the file {file} from python. This is usually an issue with the specified file path"
                )

            self.avl.avl()
            self.avl.loadgeo(geo_file)

            if  mass_file is not None:
                self.avl.loadmass(mass_file)

        else:
            raise ValueError("neither a geometry file or aircraft object was given")

        return

    def addConstraint(self, variable, val):
        self.__exe = False

        options = {
            "alpha": ["A", "A "],
            "beta": ["B", "B "],
            "roll rate": ["R", "R "],
            "pitch rate": ["P", "P "],
            "yaw rate": ["Y", "Y"],
            "D1": ["D1", "PM "],
            "D2": ["D2", "YM "],
            "D3": ["D3", "RM "],
        }

        if not (variable in options):
            raise ValueError(
                f"constraint variable `{variable}` not a valid option. Must be one of the following {[key for key in options]}."
            )

        self.avl.conset(options[variable][0], (options[variable][1] + str(val) + " \n"))
        return

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

        return

    def executeRun(self):
        self.__exe = True

        self.avl.oper()

        self.alpha = np.append(self.alpha, float(self.avl.case_r.alfa))  # *(180.0/np.pi) # returend in radians)
        self.CL = np.append(self.CL, float(self.avl.case_r.cltot))
        self.CD = np.append(self.CD, float(self.avl.case_r.cdtot))  # = np append(self.avl.case_r.cdvtot)  for total viscous)
        self.CDV = np.append(self.CDV, float(self.avl.case_r.cdvtot))  # = np append(self.avl.case_r.cdvtot)  for total viscous)
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
        sec_Chords = []
        sec_Yles = []
        end = 0
        for i in range(0, len(num_strips)):
            start = end
            end = start + num_strips[i]
            sec_CLs.append(self.avl.strp_r.clstrp[start:end])
            sec_CDs.append(self.avl.strp_r.cdstrp[start:end])
            sec_Chords.append(self.avl.strp_r.chord[start:end])
            sec_Yles.append(self.avl.strp_r.rle[1][start:end])

        self.sec_CL.append(sec_CLs)
        self.sec_CD.append(sec_CDs)
        self.sec_Chord.append(sec_Chords)
        self.sec_Yle.append(sec_Yles)

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

        

    def resetData(self):
        self.__exe = False

        self.alpha = np.zeros(0)
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
        self.sec_Chord = []
        self.sec_Yle = []

        self.surf_CL = np.empty(0)
        self.surf_CD = np.empty(0)

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
