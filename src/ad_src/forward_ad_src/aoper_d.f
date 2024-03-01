C        Generated by TAPENADE     (INRIA, Ecuador team)
C  Tapenade 3.16 (develop) - 15 Jan 2021 14:26
C
C  Differentiation of calcst in forward (tangent) mode (with options i4 dr8 r8):
C   variations   of useful results: crsax cnsax
C   with respect to varying inputs: alfa crtot cntot
C GETFILE
C
C
C
C ==================== addition wrapper helper programs ==============
C
      SUBROUTINE CALCST_D()
      INCLUDE 'AVL.INC'
      INCLUDE 'AVL_ad_seeds.inc'
      CHARACTER*50 satype
      REAL wrot_rx(3), wrot_rz(3), wrot_a(3)
      REAL crsax_u(numax), cmsax_u(numax), cnsax_u(numax), crsax_d(ndmax
     +     ), cmsax_d(ndmax), cnsax_d(ndmax), crsax_g(ngmax), cmsax_g(
     +     ngmax), cnsax_g(ngmax)
      REAL dir
      EXTERNAL GETSA
      REAL ca
      REAL ca_diff
      INTRINSIC COS
      REAL sa
      REAL sa_diff
      INTRINSIC SIN
      REAL rx
      REAL ry
      REAL rz
      REAL crsax_a
      REAL cnsax_a
      INTEGER k
      REAL cl_al
      REAL cm_al
      REAL cr_rz
      REAL cn_be
      INTRINSIC ABS
      REAL bb
      REAL cr_be
      REAL cn_rz
      REAL abs0
C
      CALL GETSA(lnasa_sa, satype, dir)
C CALL VINFAB
C CALL AERO
C
C---- set freestream velocity components from alpha, beta
C
C---- calculate forces and sensitivities
C
C---- set stability-axes rates (RX,RY,RZ) in terms of body-axes rates
      ca_diff = -(SIN(alfa)*alfa_diff)
      ca = COS(alfa)
      sa_diff = COS(alfa)*alfa_diff
      sa = SIN(alfa)
C
      rx = (wrot(1)*ca+wrot(3)*sa)*dir
      ry = wrot(2)
      rz = (wrot(3)*ca-wrot(1)*sa)*dir
C
C---- now vice-versa, and set sensitivities (which is what's really needed)
Cc    WROT(1)    =  RX*CA - RZ*SA
Cc    WROT(2)    =  RY
Cc    WROT(3)    =  RZ*CA + RX*SA
C
      wrot_rx(1) = ca*dir
      wrot_rx(2) = 0.
      wrot_rx(3) = sa*dir
C
      wrot_rz(1) = -(sa*dir)
      wrot_rz(2) = 0.
      wrot_rz(3) = ca*dir
C
C!! = -WROT(3)
      wrot_a(1) = -(rx*sa) - rz*ca
      wrot_a(2) = 0.
C!! =  WROT(1)
      wrot_a(3) = -(rz*sa) + rx*ca
C
C
      crsax_diff = dir*(ca*crtot_diff+crtot*ca_diff+sa*cntot_diff+cntot*
     +  sa_diff)
      crsax = dir*(crtot*ca+cntot*sa)
      cmsax = cmtot
      cnsax_diff = dir*(ca*cntot_diff+cntot*ca_diff-sa*crtot_diff-crtot*
     +  sa_diff)
      cnsax = dir*(cntot*ca-crtot*sa)
      crsax_a = -(crtot*sa) + cntot*ca
      cnsax_a = -(cntot*sa) - crtot*ca
C
      DO k=1,6
        crsax_u(k) = crtot_u(k)*ca + cntot_u(k)*sa
        cmsax_u(k) = cmtot_u(k)
        cnsax_u(k) = cntot_u(k)*ca - crtot_u(k)*sa
      ENDDO
C
      DO k=1,ncontrol
        crsax_d(k) = crtot_d(k)*ca + cntot_d(k)*sa
        cmsax_d(k) = cmtot_d(k)
        cnsax_d(k) = cntot_d(k)*ca - crtot_d(k)*sa
      ENDDO
C
      DO k=1,ndesign
        crsax_g(k) = crtot_g(k)*ca + cntot_g(k)*sa
        cmsax_g(k) = cmtot_g(k)
        cnsax_g(k) = cntot_g(k)*ca - crtot_g(k)*sa
      ENDDO
C
C
C---- set force derivatives in stability axes
      cltot_al = cltot_u(1)*vinf_a(1) + cltot_u(4)*wrot_a(1) + cltot_u(2
     +  )*vinf_a(2) + cltot_u(5)*wrot_a(2) + cltot_u(3)*vinf_a(3) + 
     +  cltot_u(6)*wrot_a(3) + cltot_a
      cltot_be = cltot_u(1)*vinf_b(1) + cltot_u(2)*vinf_b(2) + cltot_u(3
     +  )*vinf_b(3)
      cltot_rx = cltot_u(4)*wrot_rx(1) + cltot_u(6)*wrot_rx(3)
      cltot_ry = cltot_u(5)
      cltot_rz = cltot_u(6)*wrot_rz(3) + cltot_u(4)*wrot_rz(1)
C
      cdtot_al = cdtot_u(1)*vinf_a(1) + cdtot_u(4)*wrot_a(1) + cdtot_u(2
     +  )*vinf_a(2) + cdtot_u(5)*wrot_a(2) + cdtot_u(3)*vinf_a(3) + 
     +  cdtot_u(6)*wrot_a(3) + cdtot_a
      cdtot_be = cdtot_u(1)*vinf_b(1) + cdtot_u(2)*vinf_b(2) + cdtot_u(3
     +  )*vinf_b(3)
      cdtot_rx = cdtot_u(4)*wrot_rx(1) + cdtot_u(6)*wrot_rx(3)
      cdtot_ry = cdtot_u(5)
      cdtot_rz = cdtot_u(6)*wrot_rz(3) + cdtot_u(4)*wrot_rz(1)
C
      cytot_al = cytot_u(1)*vinf_a(1) + cytot_u(4)*wrot_a(1) + cytot_u(2
     +  )*vinf_a(2) + cytot_u(5)*wrot_a(2) + cytot_u(3)*vinf_a(3) + 
     +  cytot_u(6)*wrot_a(3)
      cytot_be = cytot_u(1)*vinf_b(1) + cytot_u(2)*vinf_b(2) + cytot_u(3
     +  )*vinf_b(3)
      cytot_rx = cytot_u(4)*wrot_rx(1) + cytot_u(6)*wrot_rx(3)
      cytot_ry = cytot_u(5)
      cytot_rz = cytot_u(6)*wrot_rz(3) + cytot_u(4)*wrot_rz(1)
C
      crtot_al = crsax_u(1)*vinf_a(1) + crsax_u(4)*wrot_a(1) + crsax_u(2
     +  )*vinf_a(2) + crsax_u(5)*wrot_a(2) + crsax_u(3)*vinf_a(3) + 
     +  crsax_u(6)*wrot_a(3) + crsax_a
      crtot_be = crsax_u(1)*vinf_b(1) + crsax_u(2)*vinf_b(2) + crsax_u(3
     +  )*vinf_b(3)
      crtot_rx = crsax_u(4)*wrot_rx(1) + crsax_u(6)*wrot_rx(3)
      crtot_ry = crsax_u(5)
      crtot_rz = crsax_u(6)*wrot_rz(3) + crsax_u(4)*wrot_rz(1)
C
      cmtot_al = cmsax_u(1)*vinf_a(1) + cmsax_u(4)*wrot_a(1) + cmsax_u(2
     +  )*vinf_a(2) + cmsax_u(5)*wrot_a(2) + cmsax_u(3)*vinf_a(3) + 
     +  cmsax_u(6)*wrot_a(3)
      cmtot_be = cmsax_u(1)*vinf_b(1) + cmsax_u(2)*vinf_b(2) + cmsax_u(3
     +  )*vinf_b(3)
      cmtot_rx = cmsax_u(4)*wrot_rx(1) + cmsax_u(6)*wrot_rx(3)
      cmtot_ry = cmsax_u(5)
      cmtot_rz = cmsax_u(6)*wrot_rz(3) + cmsax_u(4)*wrot_rz(1)
C
      cntot_al = cnsax_u(1)*vinf_a(1) + cnsax_u(4)*wrot_a(1) + cnsax_u(2
     +  )*vinf_a(2) + cnsax_u(5)*wrot_a(2) + cnsax_u(3)*vinf_a(3) + 
     +  cnsax_u(6)*wrot_a(3) + cnsax_a
      cntot_be = cnsax_u(1)*vinf_b(1) + cnsax_u(2)*vinf_b(2) + cnsax_u(3
     +  )*vinf_b(3)
      cntot_rx = cnsax_u(4)*wrot_rx(1) + cnsax_u(6)*wrot_rx(3)
      cntot_ry = cnsax_u(5)
      cntot_rz = cnsax_u(6)*wrot_rz(3) + cnsax_u(4)*wrot_rz(1)
C
C        
C
      IF (cl_al .NE. 0.0) xnp = xyzref(1) - cref*cm_al/cl_al
      IF (cr_rz*cn_be .GE. 0.) THEN
        abs0 = cr_rz*cn_be
      ELSE
        abs0 = -(cr_rz*cn_be)
      END IF
C apply the facotors to the outputs as done in the print statements of DERMATS
C
      IF (abs0 .GT. 0.0001) THEN
        bb = cr_be*cn_rz/(cr_rz*cn_be)
C        WRITE(LU,8402) BB 
C  8402  FORMAT(/' Clb Cnr / Clr Cnb  =', F11.6,
C      &    '    (  > 1 if spirally stable )')
      END IF
      crtot_al = dir*crtot_al
      crtot_be = dir*crtot_be
      cntot_al = dir*cntot_al
      cntot_be = dir*cntot_be
      cltot_rx = cltot_rx*2.0/bref
      cltot_ry = cltot_ry*2.0/cref
      cltot_rz = cltot_rz*2.0/bref
      cdtot_rx = cdtot_rx*2.0/bref
      cdtot_ry = cdtot_ry*2.0/cref
      cdtot_rz = cdtot_rz*2.0/bref
      cytot_rx = cytot_rx*2.0/bref
      cytot_ry = cytot_ry*2.0/cref
      cytot_rz = cytot_rz*2.0/bref
      crtot_rx = dir*crtot_rx*2.0/bref
      crtot_ry = dir*crtot_ry*2.0/cref
      crtot_rz = dir*crtot_rz*2.0/bref
      cmtot_rx = cmtot_rx*2.0/bref
      cmtot_ry = cmtot_ry*2.0/cref
      cmtot_rz = cmtot_rz*2.0/bref
      cntot_rx = dir*cntot_rx*2.0/bref
      cntot_ry = dir*cntot_ry*2.0/cref
      cntot_rz = dir*cntot_rz*2.0/bref
C
C C
      RETURN
C
C        WRITE(*,8401) XNP
 8401 FORMAT(/'Neutral point  Xnp =',f11.6)
      END

C  Differentiation of get_res in forward (tangent) mode (with options i4 dr8 r8):
C   variations   of useful results: alfa beta vinf delcon xyzref
C                mach cdref res res_d wv_gam
C   with respect to varying inputs: ysym zsym parval conval rv1
C                rv2 rv rc chordv enc enc_d gam gam_d
C   RW status of diff variables: ysym:in zsym:in alfa:out beta:out
C                vinf:out parval:in conval:in delcon:out xyzref:out
C                mach:out cdref:out rv1:in rv2:in rv:in rc:in chordv:in
C                enc:in enc_d:in gam:in gam_d:in res:out res_d:out
C                wv_gam:out
C
C ======================== res and Adjoint for GAM ========      
      SUBROUTINE GET_RES_D()
C
      INCLUDE 'AVL.INC'
      INCLUDE 'AVL_ad_seeds.inc'
      INTEGER i, ic
      REAL rhs_d(nvor)
      REAL rhs_d_diff(nvor)
      REAL betm
      REAL betm_diff
      INTRINSIC SQRT
      REAL(kind=8) arg1
      REAL(kind=8) arg1_diff
      REAL(kind=avl_real) temp
      INTEGER ii1
      INTEGER ii2
      CALL SET_PAR_AND_CONS_D(nitmax, irun)
C Do not use this routine in the sovler
C IF(.NOT.LAIC) THEN
C      CALL build_AIC
C end if
C---  
      CALL BUILD_AIC_D()
      amach_diff = mach_diff
      amach = mach
      arg1_diff = -(2*amach*amach_diff)
      arg1 = 1.0 - amach**2
      temp = SQRT(arg1)
      IF (arg1 .EQ. 0.D0) THEN
        betm_diff = 0.D0
      ELSE
        betm_diff = arg1_diff/(2.0*temp)
      END IF
      betm = temp
      CALL VVOR_D(betm, betm_diff, iysym, ysym, ysym_diff, izsym, zsym, 
     +            zsym_diff, vrcore, nvor, rv1, rv1_diff, rv2, rv2_diff
     +            , nsurfv, chordv, chordv_diff, nvor, rv, rv_diff, 
     +            nsurfv, .true., wv_gam, wv_gam_diff, nvmax)
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
      DO ii1=1,ndmax
        DO ii2=1,nvmax
          res_d_diff(ii2, ii1) = 0.D0
        ENDDO
      ENDDO
      DO ii1=1,nvor
        rhs_d_diff(ii1) = 0.D0
      ENDDO
C---- Setup variational BC's at the control points
      DO ic=1,ncontrol
C------ don't bother if this control variable is undefined
        IF (lcondef(ic)) THEN
          CALL MAT_PROD_D(aicn, aicn_diff, gam_d(:, ic), gam_d_diff(:, 
     +                    ic), nvor, res_d(:, ic), res_d_diff(:, ic))
C  RHS_D(:) = 0.D0
          CALL SET_GAM_D_RHS_D(ic, enc_d, enc_d_diff, rhs_d, rhs_d_diff)
          DO i=1,nvor
            res_d_diff(i, ic) = res_d_diff(i, ic) - rhs_d_diff(i)
            res_d(i, ic) = res_d(i, ic) - rhs_d(i)
          ENDDO
        END IF
      ENDDO
      END

