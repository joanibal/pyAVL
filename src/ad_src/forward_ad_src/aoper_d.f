C        Generated by TAPENADE     (INRIA, Ecuador team)
C  Tapenade 3.16 (develop) - 15 Jan 2021 14:26
C
C  Differentiation of get_res in forward (tangent) mode (with options i4 dr8 r8):
C   variations   of useful results: alfa beta vinf xyzref mach
C                res wv_gam
C   with respect to varying inputs: ysym zsym conval xyzref rv1
C                rv2 rv rc chordv enc gam
C   RW status of diff variables: ysym:in zsym:in alfa:out beta:out
C                vinf:out conval:in xyzref:in-out mach:out rv1:in
C                rv2:in rv:in rc:in chordv:in enc:in gam:in res:out
C                wv_gam:out
Csubroutine solve_rhs
      SUBROUTINE GET_RES_D()
      INCLUDE 'AVL.INC'
      INCLUDE 'AVL_ad_seeds.inc'
      REAL betm
      INTRINSIC SQRT
      INTEGER i
      REAL(kind=8) arg1
      CALL SET_PAR_AND_CONS_D(nitmax, irun)
C Do not use this routine in the sovler
C IF(.NOT.LAIC) THEN
C      CALL build_AIC
C end if
C CALL SETUP
C---  
      CALL BUILD_AIC_D()
      amach = mach
      arg1 = 1.0 - amach**2
      betm = SQRT(arg1)
      CALL VVOR_D(betm, iysym, ysym, ysym_diff, izsym, zsym, zsym_diff, 
     +            vrcore, nvor, rv1, rv1_diff, rv2, rv2_diff, nsurfv, 
     +            chordv, chordv_diff, nvor, rv, rv_diff, nsurfv, .true.
     +            , wv_gam, wv_gam_diff, nvmax)
C---- set VINF() vector from initial ALFA,BETA
      CALL VINFAB_D()
      CALL SET_VEL_RHS_D()
C
      CALL MAT_PROD_D(aicn, aicn_diff, gam, gam_diff, nvor, res, 
     +                res_diff)
C---- add the RHS vector to the residual
      DO i=1,nvor
        res_diff(i) = res_diff(i) - rhs_diff(i)
        res(i) = res(i) - rhs(i)
      ENDDO
      mach_diff = 0.D0
      END

