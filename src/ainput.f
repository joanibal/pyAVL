C***********************************************************************
C    Module:  ainput.f
C 
C    Copyright (C) 2002 Mark Drela, Harold Youngren
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************

      SUBROUTINE INPUT(LUN,FNAME,FERR)
C---------------------------------------------------------
C     Reads an processes an AVL configuration input file
C---------------------------------------------------------
      INCLUDE 'AVL.INC'
C
      CHARACTER*(*) FNAME
      LOGICAL FERR
C
      CHARACTER*4  KEYWD
      CHARACTER*120 CNAME, ANAME
      CHARACTER*128 LINE
      LOGICAL  LHINGE
C
      REAL CLX(3), CDX(3)
C
      PARAMETER (NWRK=NSMAX)
      REAL BRADY(NWRK), BRADZ(NWRK)
C
      REAL XB(IBX), YB(IBX)
      REAL XIN(IBX), YIN(IBX), TIN(IBX)
      REAL XBOD(IBX), YBOD(IBX), TBOD(IBX)

C
C---- max number of control or design variable declaration lines per section
C
C
C
C
      REAL    RINPUT(10)
      INTEGER IINPUT(10)
      LOGICAL ERROR
C
C----------------------------------------------------
      FERR = .FALSE.
C       
      OPEN(UNIT=LUN,FILE=FNAME,STATUS='OLD',ERR=3)
      GO TO 6
C
 3    CONTINUE
      CALL STRIP(FNAME,NFN)
      WRITE(*,*) 
      WRITE(*,*) '** Open error on file: ', FNAME(1:NFN)
      FNAME = FNAME(1:NFN) // '.avl'
      WRITE(*,*) '   Trying alternative: ', FNAME(1:NFN+4)
      OPEN(UNIT=LUN,FILE=FNAME,STATUS='OLD',ERR=4)
      GO TO 6
C
 4    CONTINUE
      CALL STRIP(FNAME,NFN)
      WRITE(*,*) 
      WRITE(*,*) '** Open error on file: ', FNAME(1:NFN)
      FERR = .TRUE.
      RETURN
C
C----------------------------------------------------
 6    CONTINUE
      CALL STRIP(FNAME,NFN)
      if(LVERBOSE) then 
        WRITE(*,*) 
        WRITE(*,*) 'Reading file: ', FNAME(1:NFN), '  ...'
      end if
C
      DCL_A0 = 0.
      DCM_A0 = 0.
      DCL_U0 = 0.
      DCM_U0 = 0.
C
C---- initialize all entity counters
      NSEC = 0
C
      NSURF = 0
      NVOR = 0
      NSTRIP = 0
C
      NBODY = 0
      NLNODE = 0
C
      NCONTROL = 0
      NDESIGN = 0
C
      LVISC  = .FALSE.
C
C---- initialize counters and active-entity indicators
      ISURF = 0
      IBODY = 0
C
C---- initialize input-file line counter
      ILINE = 0
C
C------------------------------------------------------------------------------
C---- start reading file
C
C---------------------------------------------------
      CALL RDLINE(LUN,LINE,NLINE,ILINE)
      TITLE = LINE(1:NLINE)
C
      if(LVERBOSE) then 
        WRITE(*,*) 
        WRITE(*,1001) TITLE(1:60)
      end if
 1001 FORMAT(/' Configuration: ', A)
C
C
C---------------------------------------------------
      CALL RDLINE(LUN,LINE,NLINE,ILINE)
      READ (LINE,*,ERR=990) MACH0
C
C---------------------------------------------------
      CALL RDLINE(LUN,LINE,NLINE,ILINE)
      READ (LINE,*,ERR=990) IYSYM, IZSYM, ZSYM
C
C---- y-symmetry plane hard-wired at y=0
      YSYM = 0.
C
      IF(IYSYM.GT.0) IYSYM =  1
      IF(IYSYM.LT.0) IYSYM = -1
      IF(IZSYM.GT.0) IZSYM =  1
      IF(IZSYM.LT.0) IZSYM = -1
C
C---------------------------------------------------
      CALL RDLINE(LUN,LINE,NLINE,ILINE)
      READ (LINE,*,ERR=990) SREF,CREF,BREF
C
      IF(SREF .LE. 0.) SREF = 1.
      IF(CREF .LE. 0.) CREF = 1.
      IF(BREF .LE. 0.) BREF = 1.
C
C---------------------------------------------------
      CALL RDLINE(LUN,LINE,NLINE,ILINE)
      READ (LINE,*,ERR=990) XYZREF0(1), XYZREF0(2), XYZREF0(3)
C
C---------------------------------------------------
C---- try to read CD data which may or may not be present
      CALL RDLINE(LUN,LINE,NLINE,ILINE)
      READ (LINE,*,ERR=8) CDREF0
C
C---- drag data was read OK... just keep going
      GO TO 10
C
 8    CONTINUE
C---- read error occurred (drag data wasn't present)...
C      ... interpret the line as keyword
      CDREF0 = 0.
      GO TO 11
C
C==============================================================================
C---- start of keyword-interpretation loop
 10   CONTINUE
      CALL RDLINE(LUN,LINE,NLINE,ILINE)
C
 11   CONTINUE
      KEYWD = LINE(1:4)
      CALL TOUPER(KEYWD)
C
C===========================================================================
      IF    (KEYWD.EQ.'EOF ') THEN
C------ end of file... clean up loose ends
C
        IF(ISURF.NE.0) THEN
C------- "old" surface is still active, so build it before finishing
         CALL MAKESURF(ISURF)
C
         IF(LDUPL(ISURF)) THEN
          CALL SDUPL(ISURF,YDUPL(ISURF),'YDUP')
          NSURF = NSURF + 1
         ENDIF
C
         ISURF = 0
        ENDIF
C
        IF(IBODY.NE.0) THEN
C------- "old" body is still active, so build it before finishing
         CALL MAKEBODY(IBODY, IBX,
     &       NVB, BSPACE,
     &       XBOD,YBOD,TBOD,NBOD)
C
         IF(LDUPL_B(IBODY)) THEN
          CALL BDUPL(IBODY,YDUPL_B(IBODY),'YDUP')
         ENDIF
C
         IBODY = 0
        ENDIF
C
C------ go finish up
        GO TO 900
C
C===========================================================================
      ELSEIF(KEYWD.EQ.'SURF') THEN
C------ new surface is about to start
C
        IF(ISURF.NE.0) THEN
C------- "old" surface is still active, so build it before starting new one
         CALL MAKESURF(ISURF)
C
         IF(LDUPL(ISURF)) THEN
          CALL SDUPL(ISURF,YDUPL(ISURF),'YDUP')
          NSURF = NSURF + 1
         ENDIF
C
         ISURF = 0
        ENDIF
C
        IF(IBODY.NE.0) THEN
C------- "old" body is still active, so build it before finishing
         CALL MAKEBODY(IBODY, IBX,
     &       NVB, BSPACE,
     &       XBOD,YBOD,TBOD,NBOD)
C
         IF(LDUPL_B(IBODY)) THEN
          CALL BDUPL(IBODY,YDUPL_B(IBODY),'YDUP')
         ENDIF
C
         IBODY = 0
        ENDIF
C
C------ new surface  (ISURF.ne.0 denotes surface accumulation is active)
        NSURF = NSURF + 1
        ISURF = MIN( NSURF , NSMAX )
C
C------ default surface component index is just the surface number
        LSCOMP(ISURF) = ISURF
C
C------ clear indices for accumulation
        NSEC(ISURF) = 0
        ISEC = 0
C
C------ set surface defaults
        LDUPL(ISURF)  = .FALSE.   
        LHINGE = .FALSE.

C------ assume this will be a conventional loaded surface
        LFWAKE(ISURF) = .TRUE.
        LFALBE(ISURF) = .TRUE.
        LFLOAD(ISURF) = .TRUE.

        XYZSCAL(1, ISURF) = 1.0
        XYZSCAL(2, ISURF) = 1.0
        XYZSCAL(3, ISURF) = 1.0
        XYZTRAN(1, ISURF) = 0.
        XYZTRAN(2, ISURF) = 0.
        XYZTRAN(3, ISURF) = 0.
        ADDINC(ISURF) = 0.
C
        NCVAR = 0
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
        STITLE(ISURF) = LINE(1:NLINE)
        
        if(lverbose)then
                WRITE(*,*)
                WRITE(*,*) '  Building surface: ', STITLE(ISURF)
        end if 
                
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
        NINPUT = 4
        CALL GETFLT(LINE,RINPUT,NINPUT,ERROR)
        IF(ERROR .OR. NINPUT.LT.2) GO TO 990
C
        NVC(ISURF) = INT( RINPUT(1) + 0.001 )
        CSPACE(ISURF) = RINPUT(2)
C
        IF(NINPUT.GE.4) THEN
         NVS(ISURF) = INT( RINPUT(3) + 0.001 )
         SSPACE(ISURF) = RINPUT(4)
        ELSE
         NVS(ISURF) = 0
         SSPACE(ISURF) = 0.0
        ENDIF
C
C===========================================================================
      ELSEIF(KEYWD.EQ.'BODY') THEN
C------ new body is about to start
C
        IF(ISURF.NE.0) THEN
C------- "old" surface is still active, so build it before starting new one     
         CALL MAKESURF(ISURF)
C
         IF(LDUPL(ISURF)) THEN
          CALL SDUPL(ISURF,YDUPL(ISURF),'YDUP')
          NSURF = NSURF + 1
         ENDIF
C
         ISURF = 0
        ENDIF
C
        IF(IBODY.NE.0) THEN
C------- "old" body is still active, so build it before finishing
         CALL MAKEBODY(IBODY, IBX,
     &       NVB, BSPACE,
     &       XBOD,YBOD,TBOD,NBOD)
C
         IF(LDUPL_B(IBODY)) THEN
          CALL BDUPL(IBODY,YDUPL_B(IBODY),'YDUP')
         ENDIF
C
         IBODY = 0
        ENDIF
C
C------ new body  (IBODY.ne.0 denotes body accumulation is active)
        NBODY = NBODY + 1
        IBODY = MIN( NBODY , NBMAX )
C
        NSEC = 0
        ISEC = 0
        NBOD = 0
C
        NIN = 0
C
        LDUPL_B(IBODY)  = .FALSE.
C       TODO: add XYZ vecotors for bodys
        XYZSCAL_B(1, IBODY) = 1.0
        XYZSCAL_B(2, IBODY) = 1.0
        XYZSCAL_B(3, IBODY) = 1.0
        XYZTRAN_B(1, IBODY) = 0.
        XYZTRAN_B(2, IBODY) = 0.
        XYZTRAN_B(3, IBODY) = 0.
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
        BTITLE(IBODY) = LINE(1:NLINE)
C         WRITE(*,*)
C         WRITE(*,*) '  Building body: ', BTITLE(IBODY)
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
        READ(LINE,*,ERR=990) NVB, BSPACE
C
C===========================================================================
      ELSEIF(KEYWD.EQ.'YDUP') THEN
C------ this surface is to be duplicated with an image surface
        IF    (ISURF.NE.0) THEN
CC         WRITE(*,*) '  + duplicate surface ',STITLE(ISURF)
        ELSEIF(IBODY.NE.0) THEN
CC         WRITE(*,*) '  + duplicate body ',BTITLE(IBODY)
        ELSE
         WRITE(*,9000) '** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)    '** No active surface or body to duplicate'
         GO TO 10
        ENDIF
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
        
C
        IF    (ISURF.NE.0) THEN
CC         WRITE(*,*) '  + duplicate surface ',STITLE(ISURF)
           LDUPL(ISURF) = .TRUE.
           READ (LINE,*,ERR=990) YDUPL(ISURF)
        ELSEIF(IBODY.NE.0) THEN
CC         WRITE(*,*) '  + duplicate body ',BTITLE(IBODY)
           LDUPL_B(IBODY) = .TRUE.
           READ (LINE,*,ERR=990) YDUPL_B(IBODY)
        ENDIF
C
        IF(IYSYM.NE.0) THEN
         WRITE(*,*)'** Warning: Y-duplicate AND Y-sym specified'
        ENDIF
C
C===========================================================================
      ELSEIF (KEYWD.EQ.'INDE' .OR. KEYWD.EQ.'COMP') THEN
C------ set component index for surface (may be lumped together)
        IF(ISURF.EQ.0) THEN
         WRITE(*,9000) '** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)    '** No active surface for component index'
         GO TO 10
        ENDIF
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
        READ (LINE,*,ERR=990) LSCOMP(ISURF)
C
C===========================================================================
      ELSEIF (KEYWD.EQ.'SCAL') THEN
C------ read scaling factors
        IF(ISURF.EQ.0 .AND. IBODY.EQ.0) THEN
         WRITE(*,9000) '** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)    '** No active surface or body for scaling'
         GO TO 10
        ENDIF
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
        
        IF    (ISURF.NE.0) THEN
                READ (LINE,*,ERR=990) XYZSCAL(1, ISURF),
     &                                XYZSCAL(2, ISURF),
     &                                XYZSCAL(3, ISURF)

        ELSEIF(IBODY.NE.0) THEN
                READ (LINE,*,ERR=990) XYZSCAL_B(1, IBODY),
     &                                XYZSCAL_B(2, IBODY),
     &                                XYZSCAL_B(3, IBODY)
        END IF 
                                
C
C===========================================================================
      ELSEIF (KEYWD.EQ.'TRAN') THEN
C------ read translation vector
        IF(ISURF.EQ.0 .AND. IBODY.EQ.0) THEN
         WRITE(*,9000) '** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)    '** No active surface or body for translation'
         GO TO 10
        ENDIF
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)

        IF    (ISURF.NE.0) THEN
               READ (LINE,*,ERR=990) XYZTRAN(1, ISURF),
     &                               XYZTRAN(2, ISURF),
     &                               XYZTRAN(3, ISURF)

        ELSEIF(IBODY.NE.0) THEN
                READ (LINE,*,ERR=990) XYZTRAN_B(1, IBODY),
     &                                XYZTRAN_B(2, IBODY),
     &                                XYZTRAN_B(3, IBODY)
        END IF 

C===========================================================================
      ELSEIF (KEYWD.EQ.'ANGL') THEN
C------ read surface angle change
        IF(ISURF.EQ.0 .AND. IBODY.EQ.0) THEN
         WRITE(*,9000) '** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)    '** No active surface or body for rotation'
         GO TO 10
        ENDIF
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
        READ (LINE,*,ERR=990) ADDINC(ISURF)
        ! ADDINC(ISURF) = ADDINC(ISURF)*DTR
        
C
C===========================================================================
      ELSEIF (KEYWD.EQ.'NOWA') THEN
C------ disable wake shedding for this surface
        IF(ISURF.EQ.0) THEN
         WRITE(*,9000)'** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)   '** No active surface for wake-shedding flag'
         GO TO 10
        ENDIF
C
        LFWAKE(ISURF) = .FALSE.

C===========================================================================
      ELSEIF (KEYWD.EQ.'NOAL') THEN
C------ disable freestream angles for this surface
        IF(ISURF.EQ.0) THEN
         WRITE(*,9000)'** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)   '** No active surface for freestream-angles flag'
         GO TO 10
        ENDIF
C
        LFALBE(ISURF) = .FALSE.

C===========================================================================
      ELSEIF (KEYWD.EQ.'NOLO') THEN
C------ disable total-load contributions for this surface
        IF(ISURF.EQ.0) THEN
         WRITE(*,9000)'** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)   '** No active surface for load-disable flag'
         GO TO 10
        ENDIF
C
        LFLOAD(ISURF) = .FALSE.

C===========================================================================
      ELSEIF (KEYWD.EQ.'BSEC') THEN
C------ read body section
C
        IF(IBODY.EQ.0) THEN
         WRITE(*,9000) '** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)    '** No active body for this section'
         GO TO 10
        ENDIF
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
C
C------ store section data for current body
        NSEC_B(IBODY) = NSEC_B(IBODY) + 1
        ISEC = MIN( NSEC_B(IBODY) , NWRK )
C
        NINPUT = 5
        CALL GETFLT(LINE,RINPUT,NINPUT,ERROR)
        IF(ERROR .OR. NINPUT.LT.4) GO TO 990
C
        XYZLES_B(1,ISEC,IBODY) = RINPUT(1)
        XYZLES_B(2,ISEC,IBODY) = RINPUT(2)
        XYZLES_B(3,ISEC,IBODY) = RINPUT(3)
        BRADY(ISEC) = RINPUT(4)
C
        IF(NINPUT.GE.5) THEN
         BRADZ(ISEC) = RINPUT(5)
        ELSE
         BRADZ(ISEC) = BRADY(ISEC)
        ENDIF
C
C===========================================================================
      ELSEIF (KEYWD.EQ.'SECT') THEN
C------ read surface section
C
        IF(ISURF.EQ.0) THEN
         WRITE(*,9000) '** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)    '** No active surface for this section'
         GO TO 10
        ENDIF
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
C
C------ store section data for current surface
        NSEC(ISURF) = NSEC(ISURF) + 1
        ISEC = MIN( NSEC(ISURF) , NWRK )
C
        NINPUT = 7
        CALL GETFLT(LINE,RINPUT,NINPUT,ERROR)
        IF(ERROR .OR. NINPUT.LT.5) GO TO 990
C
        XYZLES(1,ISEC, ISURF) = RINPUT(1)
        XYZLES(2,ISEC, ISURF) = RINPUT(2)
        XYZLES(3,ISEC, ISURF) = RINPUT(3)
        CHORDS(ISEC, ISURF) = RINPUT(4)
        AINCS(ISEC, ISURF)  = RINPUT(5)
C
        IF(NINPUT.GE.7) THEN
         NSPANS(ISEC, ISURF) = INT( RINPUT(6) + 0.001 )
         SSPACES(ISEC, ISURF) = RINPUT(7)
        ELSE
         NSPANS(ISEC, ISURF) = 0
         SSPACES(ISEC, ISURF) = 0.
        ENDIF
C
C------ default section...
C    ...flat camberline
        NASEC(ISEC,ISURF)   = 2
        XASEC(1,ISEC, ISURF) = 0.0
        XASEC(2,ISEC, ISURF) = 1.0
        SASEC(1,ISEC, ISURF)  = 0.0
        SASEC(2,ISEC, ISURF)  = 0.0
        TASEC(1,ISEC, ISURF)  = 0.0
        TASEC(2,ISEC, ISURF)  = 0.0
C    ...no polar data
        DO L=1, 6
          CLCDSEC(L,ISEC,ISURF) = 0.0
        END DO
C    ...no control
        NSCON(ISEC, ISURF) = 0
C    ...no design
        NSDES(ISEC, ISURF) = 0
C
C    ...unity dCL/da factor
        CLAF(ISEC,ISURF) = 1.0
C
C===========================================================================
      ELSEIF (KEYWD.EQ.'NACA') THEN 
C------ input NACA camberline
C
        IF(ISURF.EQ.0 .OR. ISEC.EQ.0) THEN
         WRITE(*,9000) '** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)    '** No active section for this airfoil'
         GO TO 10
        ENDIF
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
C
        IB = INDEX(LINE,' ')
        READ(LINE(1:IB-1),*,ERR=990) IDES
        IF(LINE(IB:NLINE).NE.' ') THEN
         READ(LINE(IB:NLINE),*,ERR=990) XFMIN, XFMAX
ccc           WRITE(*,*) '   Using data in normalized range ',XFMIN,XFMAX
        ELSE
         XFMIN = 0.
         XFMAX = 1.
        ENDIF
C
        ICAM = IDES/1000
        IPOS = (IDES-1000*ICAM)/100
        ITHK =  IDES-1000*ICAM-100*IPOS
        C = FLOAT(ICAM) / 100.0
        P = FLOAT(IPOS) / 10.0
        T = FLOAT(ITHK) / 100.0
C
        NASEC(ISEC,ISURF) = MIN( 50 , IBX )
        DO I = 1, NASEC(ISEC,ISURF)
          XF = XFMIN +
     &        (XFMAX-XFMIN)*FLOAT(I-1)/FLOAT(NASEC(ISEC,ISURF)-1)

          XASEC(I,ISEC, ISURF) = XF
          TASEC(I,ISEC, ISURF) = 
     &        (   0.29690*SQRT(XF)
     &          - 0.12600*XF
     &          - 0.35160*XF**2
     &          + 0.28430*XF**3
     &          - 0.10150*XF**4  ) * T * 10.0
          IF    (XF .LT. P) THEN
           SASEC(I,ISEC,ISURF) = C * 2.0*(P - XF) / P**2
          ELSEIF(XF .GE. P) THEN
           SASEC(I,ISEC,ISURF) = C * 2.0*(P - XF) / (1.0-P)**2
          ELSE
           SASEC(I,ISEC,ISURF) = 0.
          ENDIF
        ENDDO

        CALL NRMLIZ(NASEC(ISEC,ISURF),XASEC(1,ISEC,ISURF))
C
C===========================================================================
      ELSE IF (KEYWD.EQ.'AIRF') THEN 
C------ input y(x) for an airfoil, get camber then slopes via spline
C
        IF(ISURF.EQ.0 .OR. ISEC.EQ.0) THEN
         WRITE(*,9000) '** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)    '** No active section for this airfoil'
         GO TO 10
        ENDIF
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
C
        IB = INDEX(LINE,' ')
        READ(LINE(IB:NLINE),*,ERR=990) XFMIN, XFMAX
C
        DO I = 1, 123456
          IB = MIN(I,IBX)
C
          CALL RDLINE(LUN,LINE,NLINE,ILINE)
          NINPUT = 2
          CALL GETFLT(LINE,RINPUT,NINPUT,ERROR)
          IF(ERROR .OR. NINPUT.LT.2) THEN
           NB = IB-1
           GO TO 40
          ELSE
           XB(IB) = RINPUT(1)
           YB(IB) = RINPUT(2)
          ENDIF
        ENDDO
C
 40     CONTINUE
        IF(I.GT.IBX) THEN
         WRITE(*,*) 
     &    '*** AINPUT: Airfoil array overflow.  Increase IBX to', I
         STOP
        ENDIF
C
C------ set camber and thickness, normalized to unit chord
        NIN = MIN( 50 , IBX )
        CALL GETCAM(XB,YB,NB,XIN,YIN,TIN,NIN,.TRUE.)
C
C------ store airfoil only if surface and section are active
        NASEC(ISEC,ISURF) = NIN
        DO I = 1, NIN
          XF = XFMIN + 
     &         (XFMAX-XFMIN)*FLOAT(I-1)/FLOAT(NASEC(ISEC,ISURF)-1)
          XASEC(I,ISEC,ISURF) = XIN(1) + XF*(XIN(NIN)-XIN(1))
          CALL AKIMA(XIN,YIN,NIN,XASEC(I,ISEC,ISURF),DUMMY,
     &               SASEC(I,ISEC,ISURF))
          CALL AKIMA(XIN,TIN,NIN,XASEC(I,ISEC,ISURF),
     &               TASEC(I,ISEC,ISURF),DUMMY)
        END DO
        CALL NRMLIZ(NASEC(ISEC,ISURF),XASEC(1,ISEC,ISURF))
C
C------ go to top of keyword-reading loop, with last-read line
        GO TO 11
C
C===========================================================================
      ELSEIF (KEYWD.EQ.'AFIL') THEN 
C------ input y(x) from an airfoil coordinate file
C
        IF(ISURF.EQ.0 .OR. ISEC.EQ.0) THEN
         WRITE(*,9000) '** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)    '** No active section for this airfoil'
         GO TO 10
        ENDIF
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
C---- parse file name and optional parameters
C     double quotes checked to delimit file name to allow blanks in name 
        IDQ1 = INDEX(LINE,'"')
        IF(IDQ1.NE.0) THEN
         IDQ2 = INDEX(LINE(IDQ1+1:),'"')
         IF(IDQ2.GT.1) THEN 
           CNAME = LINE(IDQ1+1:IDQ2+IDQ1-1)
           IB = IDQ2 + IDQ1 + 1
         ELSE
           WRITE(*,9000) '** Bad quotes in file name ',ILINE,
     &                   LINE(1:NLINE)
           GO TO 10
         ENDIF
        ELSE
C---- Find blank after filename as delimiter for optional parameters
         IB = INDEX(LINE,' ')
         CNAME = LINE(1:IB)
        ENDIF
C
        IF(LINE(IB:NLINE).NE.' ') THEN
         READ(LINE(IB:NLINE),*,ERR=990) XFMIN, XFMAX
CCC         WRITE(*,*) '     Using data in normalized range ',XFMIN,XFMAX
        ELSE
         XFMIN = 0.
         XFMAX = 1.
        ENDIF
C
        CALL STRIP(CNAME,NCN)
        if(lverbose)then
        WRITE(*,*) '    Reading airfoil from file: ',CNAME(1:NCN)
        end if 
        AFILES(ISEC, ISURF) = CNAME(1:NCN)        
        NBLDS = 1
        CALL READBL(CNAME,IBX,NBLDS,XB,YB,NB,NBL,
     &               ANAME,XINL,XOUT,YBOT,YTOP)
C
        IF(NBL.EQ.0) THEN
        if(lverbose)then
        WRITE(*,*) '**   Airfoil file not found  : ',CNAME(1:NCN)
        WRITE(*,*) '**   Using default zero-camber airfoil'
        end if 

C C
         NASEC(ISEC,ISURF) = MIN( 50 , IBX )
         DO I = 1, NASEC(ISEC,ISURF)
          XASEC(I,ISEC,ISURF) = FLOAT(I-1)/FLOAT(NASEC(ISEC,ISURF)-1)
          SASEC(I,ISEC,ISURF) = 0.0
          TASEC(I,ISEC,ISURF) = 0.0
         ENDDO
C
        ELSE
C------- camber and thickness
         NIN = MIN( 50 , IBX )
         CALL GETCAM(XB,YB,NB,XIN,YIN,TIN,NIN,.TRUE.)
C
C------- camberline slopes at specified locations from spline
         NASEC(ISEC,ISURF) = NIN
         DO I = 1, NIN
           XF = XFMIN + 
     &         (XFMAX-XFMIN)*FLOAT(I-1)/FLOAT(NASEC(ISEC,ISURF)-1)
           XASEC(I,ISEC,ISURF) = XIN(1) + XF*(XIN(NIN)-XIN(1))
           CALL AKIMA(XIN,YIN,NIN,XASEC(I,ISEC,ISURF),DU
     &               MMY,SASEC(I,ISEC,ISURF))
           CALL AKIMA(XIN,TIN,NIN,XASEC(I,ISEC,ISURF),
     &               TASEC(I,ISEC,ISURF),DUMMY)
         END DO
         CALL NRMLIZ (NASEC(ISEC,ISURF),XASEC(1,ISEC,ISURF))
        ENDIF
C
C===========================================================================
      ELSEIF (KEYWD.EQ.'BFIL') THEN 
C------ input r(x) from an airfoil coordinate file
C
        IF(IBODY.EQ.0) THEN
         WRITE(*,9000) '** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)    '** No active body for this shape'
         GO TO 10
        ENDIF
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
C---- parse file name and optional parameters
C     double quotes checked to delimit file name to allow blanks in name 
        IDQ1 = INDEX(LINE,'"')
        IF(IDQ1.NE.0) THEN
         IDQ2 = INDEX(LINE(IDQ1+1:),'"')
         IF(IDQ2.GT.1) THEN 
           CNAME = LINE(IDQ1+1:IDQ2+IDQ1-1)
           IB = IDQ2 + IDQ1 + 1
         ELSE
           WRITE(*,9000) '** Bad quotes in file name ',ILINE,
     &                   LINE(1:NLINE)
           GO TO 10
         ENDIF
        ELSE
C---- Find blank after filename as delimiter for optional parameters
         IB = INDEX(LINE,' ')
         CNAME = LINE(1:IB)
        ENDIF
C
        IF(LINE(IB:NLINE).NE.' ') THEN
         READ(LINE(IB:NLINE),*,ERR=990) XFMIN, XFMAX
CCC         WRITE(*,*) '     Using data in normalized range ',XFMIN,XFMAX
        ELSE
         XFMIN = 0.
         XFMAX = 1.
        ENDIF
C
        CALL STRIP(CNAME,NCN)
        WRITE(*,*) '    Reading body shape from file: ',CNAME(1:NCN)
        NBLDS = 1
        CALL READBL(CNAME,IBX,NBLDS,XB,YB,NB,NBL,
     &               ANAME,XINL,XOUT,YBOT,YTOP)
C
C------ set thread line y, and thickness t ( = 2r)
        NBOD = MIN( 50 , IBX )
        CALL GETCAM(XB,YB,NB,XBOD,YBOD,TBOD,NBOD,.FALSE.)
C
C===========================================================================
      ELSEIF (KEYWD.EQ.'CDCL') THEN 
C------ input approximate CD(CL) polar defining data
C
        IF(ISURF.EQ.0 .OR. ISEC.EQ.0) THEN
         WRITE(*,9000) '** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)    '** No active section for this polar'
         GO TO 10
        ENDIF

        CALL RDLINE(LUN,LINE,NLINE,ILINE)
        READ(LINE,*,ERR=990) CLX(1),CDX(1),CLX(2),CDX(2),CLX(3),CDX(3)
C
        LMAX = 1
        LMIN = 1
        DO L = 2, 3 
          IF(CLX(L).GT.CLX(LMAX)) LMAX = L
          IF(CLX(L).LT.CLX(LMIN)) LMIN = L
        END DO
C
        IF(ISEC.GT.1) THEN
         IF(CLCDSEC(4,ISEC-1,ISURF).LE.0.0) THEN
          WRITE(*,*) '* AINPUT: previous section defined with no polar' 
         ENDIF
        ENDIF
C
C------ Trick: sum must be 6 so we can get the "other" index
        LMID = 6 - (LMIN+LMAX)
        CLCDSEC(1,ISEC,ISUF) = CLX(LMIN)
        CLCDSEC(2,ISEC,ISUF) = CDX(LMIN)
        CLCDSEC(3,ISEC,ISUF) = CLX(LMID)
        CLCDSEC(4,ISEC,ISUF) = CDX(LMID)
        CLCDSEC(5,ISEC,ISUF) = CLX(LMAX)
        CLCDSEC(6,ISEC,ISUF) = CDX(LMAX)
        
        if (LVERBOSE) then
                WRITE(*,1700) CLX(LMIN),CDX(LMIN),
     &                CLX(LMID),CDX(LMID),
     &                CLX(LMAX),CDX(LMAX)
1700   FORMAT('    Reading CD(CL) data for section',
     &         /'     CLneg    = ',F8.3,'  CD@CLneg = ',F10.5,
     &         /'     CL@CDmin = ',F8.3,'  CDmin    = ',F10.5,
     &         /'     CLpos    = ',F8.3,'  CD@CLpos = ',F10.5)
        ENDIF
        LVISC = .TRUE.
C
C===========================================================================
      ELSEIF (KEYWD.EQ.'CLAF') THEN 
C------ input dCL/da scaling factor
C
        IF(ISURF.EQ.0 .OR. ISEC.EQ.0) THEN
         WRITE(*,9000) '** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)    '** No active section for dCL/da factor'
         GO TO 10
        ENDIF

        CALL RDLINE(LUN,LINE,NLINE,ILINE)
        READ(LINE,*,ERR=990) CLAF(ISEC,ISURF)
C
        IF(CLAF(ISEC,ISURF) .LE. 0.0 .OR. 
     &     CLAF(ISEC,ISURF) .GE. 2.0      ) THEN
         WRITE(*,*) '** dCL/da factor must be in the range 0..2 --',
     &              ' Setting factor to 1.0'
         CLAF(ISEC,ISURF) = 1.0
        ENDIF
C
C===========================================================================
      ELSEIF (KEYWD.EQ.'CONT') THEN
C------ link section to control variables
C
        IF(ISURF.EQ.0 .OR. ISEC.EQ.0) THEN
         WRITE(*,9000) '** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)    '** No active section for this control'
         GO TO 10
        ENDIF
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
C
C
C------ increment control-declaration counter for this section
        NSCON(ISEC,ISURF) = NSCON(ISEC,ISURF) + 1
        ISCON = MIN( NSCON(ISEC,ISURF) , ICONX )
C
C------ extract control name
        NNAME = INDEX(LINE,' ') - 1
        IF(NNAME.LE.0) THEN
         WRITE(*,*) '** Bad control declaration line:  ', LINE
         STOP
        ENDIF
C
C------ see if this control variable has already been declared
        DO N = 1, NCONTROL
          IF(LINE(1:NNAME) .EQ. DNAME(N)(1:NNAME)) THEN
           ICONTROL = N
           GO TO 62
          ENDIF
        ENDDO
C
C------ new control variable... assign slot for it
        NCONTROL = NCONTROL + 1
        ICONTROL = MIN( NCONTROL , NDMAX )
        DNAME(ICONTROL) = LINE(1:NNAME)
C
 62     CONTINUE
        ICONTD(ISCON,ISEC,ISURF) = ICONTROL
C
C------ read numbers after control variable name
        NINPUT = 6
        CALL GETFLT(LINE(NNAME+1:120),RINPUT,NINPUT,ERROR)
        IF(ERROR) THEN
         WRITE(*,*) '*** Bad control data line:  ', LINE
         STOP
        ENDIF
C
        IF(NINPUT.LT.1) THEN
         GAIND(ISCON,ISEC,ISURF) = 1.0
        ELSE
         GAIND(ISCON,ISEC,ISURF) = RINPUT(1)
        ENDIF
C
        IF(NINPUT.LT.2) THEN
         XHINGED(ISCON,ISEC,ISURF) = 0.0
        ELSE
         XHINGED(ISCON,ISEC,ISURF) = RINPUT(2)
        ENDIF
C
        IF(NINPUT.LT.5) THEN
         VHINGED(1,ISCON,ISEC,ISURF) = 0.0
         VHINGED(2,ISCON,ISEC,ISURF) = 0.0
         VHINGED(3,ISCON,ISEC,ISURF) = 0.0
        ELSE
         VHINGED(1,ISCON,ISEC,ISURF) = RINPUT(3)
         VHINGED(2,ISCON,ISEC,ISURF) = RINPUT(4)
         VHINGED(3,ISCON,ISEC,ISURF) = RINPUT(5)
        ENDIF
C
        IF(NINPUT.LT.6) THEN
         REFLD(ISCON,ISEC,ISURF) = 1.0
        ELSE
         REFLD(ISCON,ISEC,ISURF) = RINPUT(6)
        ENDIF
C
C===========================================================================
      ELSEIF (KEYWD.EQ.'DESI') THEN 
C------ link section to design variable and weight
C
        IF(ISURF.EQ.0 .OR. ISEC.EQ.0) THEN
         WRITE(*,9000) '** Misplaced line', ILINE, LINE(1:NLINE)
         WRITE(*,*)    '** No active section for this design var.'
         GO TO 10
        ENDIF
C
        CALL RDLINE(LUN,LINE,NLINE,ILINE)
C
C------ increment design-declaration counter for this section
        NSDES(ISEC,ISURF) = NSDES(ISEC,ISURF) + 1
        ISDES = MIN( NSDES(ISEC,ISURF) , ICONX )
C
C------ extract design name
        NNAME = INDEX(LINE,' ') - 1
        IF(NNAME.LE.0) THEN
         WRITE(*,9000) '   *** Bad design declaration line', 
     &                  ILINE, LINE(1:NLINE)
         STOP
        ENDIF
C
C------ see if this control variable has already been declared
        DO K = 1, NDESIGN
          IF(LINE(1:NNAME) .EQ. GNAME(K)(1:NNAME)) THEN
           IDESIGN = K
           GO TO 72
          ENDIF
        ENDDO
C
        NDESIGN = NDESIGN + 1
        IDESIGN = MIN( NDESIGN , NGMAX )
        GNAME(IDESIGN) = LINE(1:NNAME)
C
 72     CONTINUE
        IDESTD(ISDES,ISEC,ISURF) = IDESIGN
C
C------ read numbers after control variable name
        NINPUT = 1
        CALL GETFLT(LINE(NNAME+1:120),RINPUT,NINPUT,ERROR)
        IF(ERROR) GO TO 990
C
        IF(NINPUT.LT.1) THEN
         GAING(ISDES,ISEC,ISURF) = 1.0
        ELSE
         GAING(ISDES,ISEC,ISURF) = RINPUT(1)
        ENDIF
C
C===========================================================================
      ELSE
C------ line not recognized or unassignable ... keep reading file
        WRITE(*,8000) ILINE, LINE(1:NLINE)
 8000   FORMAT('  * Line',I5,' ignored: ', A)
C
      ENDIF
      GO TO 10
C
C===========================================================================
C---- normal end-of-file exit point
 900  CONTINUE
      CLOSE(UNIT=LUN)
C
C*********************************************************************
C
C       WRITE (*,2018) MACH0,NBODY,NSURF,NSTRIP,NVOR
C
      IF(IYSYM.GT.0) WRITE (*,2024) YSYM
      IF(IYSYM.LT.0) WRITE (*,2025) YSYM
      IF(IZSYM.GT.0) WRITE (*,2026) ZSYM
      IF(IZSYM.LT.0) WRITE (*,2027) ZSYM
C
 2018 FORMAT (/' Mach =',F10.4,'  (default)'
     &        /' Nbody =',I4,5X,
     &         ' Nsurf =',I4,5X,
     &         ' Nstrp =',I4,5X,
     &         ' Nvor  =',I4)
 2024 FORMAT (/' Y Symmetry: Wall plane   at Ysym =',F10.4)
 2025 FORMAT (/' Y Symmetry: Free surface at Ysym =',F10.4)
 2026 FORMAT (/' Z Symmetry: Ground plane at Zsym =',F10.4)
 2027 FORMAT (/' Z Symmetry: Free surface at Zsym =',F10.4)
C
      LGEO = .TRUE.
      RETURN
C
C*********************************************************************
 990  CONTINUE
      WRITE(*,9000) '** Read error on line', ILINE, LINE(1:NLINE)
      FERR = .TRUE.
      RETURN
C
 9000 FORMAT(/ 1X,A,I5,' ...' / 1X,A)
      END ! INPUT





      SUBROUTINE RDLINE(LUN,LINE,NLINE,ILINE)
C-----------------------------------------------------------------------
C     Reads next non-comment line from logical unit LU
C     Strips off leading blanks
C     Ignores everything after and including "!"
C
C     LINE returns the line
C     NLINE returns the number of characters in non-blank portion
C
C     If e.o.f. is reached, LINE returns 'EOF'
C     If read error occurs, LINE returns 'ERR'
C-----------------------------------------------------------------------
      CHARACTER*(*) LINE
C
 1000 FORMAT(A)
   20 READ (LUN,1000,END=80,ERR=90) LINE
      ILINE = ILINE + 1
C
C---- skip comment line
      IF(INDEX('!#',LINE(1:1)) .NE. 0) GO TO 20
C
C---- skip blank line
      IF(LINE.EQ.' ') GO TO 20
C
C---- strip off leading blanks and do normal return after significant line
      CALL STRIP(LINE,NLINE)
      KEXL = INDEX(LINE(1:NLINE),'!')
      IF(KEXL.GT.1) NLINE = KEXL-1
      RETURN
C
   80 LINE = 'EOF '
      RETURN
C
   90 LINE = 'ERR '
      RETURN
      END


      SUBROUTINE TOUPER(INPUT)
      CHARACTER*(*) INPUT
C
      CHARACTER*26 LCASE, UCASE
      DATA LCASE / 'abcdefghijklmnopqrstuvwxyz' /
      DATA UCASE / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
C
      N = LEN(INPUT)
C
      DO 10 I=1, N
        K = INDEX( LCASE , INPUT(I:I) )
        IF(K.GT.0) INPUT(I:I) = UCASE(K:K)
 10   CONTINUE
C
      RETURN
      END 


