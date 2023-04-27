C        Generated by TAPENADE     (INRIA, Ecuador team)
C  Tapenade 3.16 (develop) - 15 Jan 2021 14:26
C
C  Differentiation of get_res in forward (tangent) mode (with options i4 dr8 r8):
C   variations   of useful results: res
C   with respect to varying inputs: alfa rv1 rv2 rc chordv enc
C   RW status of diff variables: alfa:in rv1:in rv2:in rc:in chordv:in
C                enc:in res:out
Csubroutine solve_rhs
      SUBROUTINE GET_RES_D()
      INCLUDE 'AVL.INC'
      INCLUDE 'AVL_ad_seeds.inc'
      INTEGER i
      INTEGER ii2
      INTEGER ii1
C---  
      IF (.NOT.laic) THEN
        CALL BUILD_AIC_D()
      ELSE
        DO ii1=1,nvmax
          DO ii2=1,nvmax
            aicn_diff(ii2, ii1) = 0.D0
          ENDDO
        ENDDO
      END IF
C---- set VINF() vector from initial ALFA,BETA
      CALL VINFAB_D()
      CALL SET_VEL_RHS_D()
C
      CALL MAT_PROD_D(aicn, aicn_diff, gam, nvor, res, res_diff)
C---- add the RHS vector to the residual
      DO i=1,nvor
        res_diff(i) = res_diff(i) - rhs_diff(i)
        res(i) = res(i) - rhs(i)
      ENDDO
      END

