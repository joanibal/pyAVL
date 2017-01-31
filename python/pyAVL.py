'''
pyAVL

pyAVL is a wrapper for Mark Drela's Xfoil code. The purpose of this
class is to provide an easy to use wrapper for avl for intergration
into other projects. 

Developers:
-----------
- Josh Anibal (JLA)

History
-------
	v. 1.0 - Initial Class Creation (JLA, 08 2016)
'''

__version__ = 1.0

# =============================================================================
# Standard Python modules
# =============================================================================

import os, sys, string, copy, pdb, time

# =============================================================================
# External Python modules
# =============================================================================

import numpy as np

# =============================================================================
# Extension modules
# =============================================================================

import pyavl as avl



class avlAnalysis():
    def __init__(self, geo_file =None, mass_file =None,  aircraft_object = None ):

        self.__exe = False

        self.alpha =[]
        self.CL = []
        self.CD = []
        self.CM = []
        self.span_eff = []

        self.elev_def = []
        self.rud_def = []

        self.velocity = []


        self.sec_CL = []
        self.sec_CD = []
        self.sec_Chord = []
        self.sec_Yle = []


    
        if not(geo_file == None):
            try:
                # check to make sure files exist 
                file = geo_file
                f = open(geo_file,'r')
                f.close()

                if not(mass_file is None):
                    file = mass_file
                    f = open(mass_file,'r')
                    f.close
            except:
                print 'ERROR:  There was an error opening the file %s'%(file)
                sys.exit(1)

            avl.avl()   
            avl.loadgeo(geo_file)
            avl.loadmass(mass_file)


        elif not(aircraft_object == None):

            pass


            # for A in [a for a in dir(obj) if not a.startswith('__')]
                # if type(plane.A) is Surface:
                    # pass

            # avl.avl()   
 
            # avl.CASE_R.SREF = aircraft.Sref
            # avl.CASE_R.CREF = aircraft.Cref
            # avl.CASE_R.BREF = aircraft.Bref


            # # avl.MASS_R.RMASS0 = aircraft.mass
            # # avl.MASS_R.XYZMASS0 = (aircraft.X_cg, aircraft.Y_cg, aircraft.Z_cg)
            # # avl.MASS_R.RINER0 = (aircraft.I)

            # # avl.MASS_R.AMASS(3,3) = aircraft.mass
            # # avl.MASS_R.AINER0(3,3) = (aircraft.I)



            # avl.CASE_L.CDREF = aircraft.CD_p

            # avl. = aircraft.Chord

            # avl. = aircraft.Xle
            # avl. = aircraft.Yle
            # avl. = aircraft.Zle

            # avl. = aircraft.Ainc

            # avl. = aircraft.Tail_Horz_Sref 
            # avl. = aircraft.Tail_Horz_Chord
            # avl. = aircraft.Tail_Horz_Bref 
            # avl. = aircraft.Tail_Horz_Ainc


            # avl. = aircraft.Tail_Vert_Sref 
            # avl. = aircraft.Tail_Vert_Chord
            # avl. = aircraft.Tail_Vert_Bref 
        else:
            print 'ERROR:  neither a geometry file or aircraft object was given'
            sys.exit(1)

        return


    def addConstraint(self,varible,val):

        self.__exe = False

        options = {'alpha':['A', 'A '],
                   'beta':['B', 'B '],
                   'roll rate':['R', 'R '],
                   'pitch rate':['P', 'P '],
                   'yaw rate':['Y', 'Y'],
                   'elevator':['D1', 'PM '],
                   'rudder': ['D2', 'YM '],
                   'aileron': ['D3', 'RM ']}


        if not(varible in options):
            print 'ERROR:  constraint varible not a valid option '
            sys.exit(1)

        avl.conset(options[varible][0],(options[varible][1] +  str(val) + ' \n'))
        return

    def addTrimCondition(self, varible, val):

        self.__exe = False

        options = {'bankAng':['B'],
                   'CL':['C'],
                   'velocity':['V'],
                   'mass':['M'],
                   'dens':['D'],
                   'G':['G'],
                   'X_cg': ['X'],
                   'Y_cg': ['Y'],
                   'Z_cg': ['Z']}

        if not(varible in options):
            print 'ERROR:  constraint varible not a valid option '
            sys.exit(1)



        avl.trmset('C1','1 ',options[varible][0],(str(val) +'  \n'))
        pass
        return


    def executeRun(self):

        self.__exe = True

        avl.oper()

        self.alpha.append(float(avl.case_r.alfa))  # *(180.0/np.pi) # returend in radians)
        self.CL.append(float(avl.case_r.cltot))
        self.CD.append(float(avl.case_r.cdtot))   # append(avl.case_r.cdvtot)  for total viscous)
        self.CM.append(float(avl.case_r.cmtot))
        self.span_eff.append(float(avl.case_r.spanef))


        self.elev_def.append(float(avl.case_r.delcon[0]))
        self.rud_def.append(float(avl.case_r.delcon[1]))

        self.velocity.append(np.asarray(avl.case_r.vinf))

        # get section properties
        NS = avl.surf_i.nj[0]

        self.sec_CL.append(np.asarray(avl.strp_r.clstrp[:NS]))
        self.sec_CD.append(np.asarray(avl.strp_r.cdstrp[:NS]))
        self.sec_Chord.append(np.asarray(avl.strp_r.chord[:NS]))
        self.sec_Yle.append(np.asarray(avl.strp_r.rle[1][:NS]))


        # print 'alfa:', avl.case_r.alfa   

        # print 'CLTOT:', avl.case_r.cltot
        # print 'CdTOT:', avl.case_r.cdtot
        # print 'CmTOT:', avl.case_r.cmtot
        # print 'Dname', avl.case_c.dname
        # print 'Delcon', avl.case_r.delcon


        

        return 

    def calcNP(self):
        # executeRun must be run first 

        if not(self.__exe):
            print 'ERROR:  executeRun most be called first'
            sys.exit(1)


        avl.calcst()
        self.NP = avl.case_r.xnp
        # print 'Xnp:', avl.case_r.xnp

        return

    

    def alphaSweep(self, start_alpha, end_alpha, increment=1):

        alphas = np.arange(start_alpha, end_alpha+increment, increment)

        for alf in alphas:
            self.addConstraint('alpha',alf)  
            self.executeRun()

        return

    def CLSweep(self, start_CL, end_CL, increment=0.1):

        CLs = np.arange(start_CL, end_CL+increment, increment)

        for cl in CLs:
            self.addTrimCondition('CL',cl)  
            self.executeRun()

        return


    def clearVals(self):
        self.__exe = False

        self.alpha =[]
        self.CL = []
        self.CD = []
        self.CM = []
        self.span_eff = []

        self.elev_def = []
        self.rud_def = []

        self.velocity = []


        self.sec_CL = []
        self.sec_CD = []
        self.sec_Chord = []
        self.sec_Yle = []


    # def calcSectionalCPDist(self):

    #     for N in [1, avl.case_i.nsurf]:

    #         J1 = avl.surf_i.jfrst(N)
    #         JN = J1 + avl.surf_i.nj - 1
    #         JINC = 1

    #         CPSCL = avl.   cpfac*  avl.case_r.cref


    #         for 


    #     return



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