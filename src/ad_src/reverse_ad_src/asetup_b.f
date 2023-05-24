C        Generated by TAPENADE     (INRIA, Ecuador team)
C  Tapenade 3.16 (develop) - 15 Jan 2021 14:26
C
C  Differentiation of build_aic in reverse (adjoint) mode (with options i4 dr8 r8):
C   gradient     of useful results: ysym zsym rv1 rv2 rc chordv
C                enc aicn
C   with respect to varying inputs: ysym zsym rv1 rv2 rc chordv
C                enc
C SETUP
C
C
      SUBROUTINE BUILD_AIC_B()
      INCLUDE 'AVL.INC'
      INCLUDE 'AVL_ad_seeds.inc'
      REAL betm
      INTRINSIC SQRT
      INTEGER j
      INTEGER i
      INTEGER n
      INTEGER j1
      INTEGER jn
      INTEGER i1
      INTEGER iv
      INTEGER jv
      INTEGER ii3
      INTEGER ii2
      INTEGER ii1
      amach = mach
      betm = SQRT(1.0 - amach**2)
      CALL VVOR(betm, iysym, ysym, izsym, zsym, vrcore, nvor, rv1, rv2, 
     +          nsurfv, chordv, nvor, rc, nsurfv, .false., wc_gam, nvmax
     +         )
C$BWD-OF II-LOOP 
      DO n=1,nsurf
        IF (.NOT.lfwake(n)) THEN
C
C------- go over TE control points on this surface
          j1 = jfrst(n)
C$AD II-LOOP
          jn = jfrst(n) + nj(n) - 1
C$BWD-OF II-LOOP 
          DO j=j1,jn
            i1 = ijfrst(j)
            iv = ijfrst(j) + nvstrp(j) - 1
            DO jv=iv,i1,-1
              aicn_diff(iv, jv) = 0.D0
            ENDDO
C$BWD-OF II-LOOP 
            DO jv=1,nvor
              aicn_diff(iv, jv) = 0.D0
            ENDDO
          ENDDO
        END IF
      ENDDO
      DO ii1=1,nvmax
        DO ii2=1,nvmax
          DO ii3=1,3
            wc_gam_diff(ii3, ii2, ii1) = 0.D0
          ENDDO
        ENDDO
      ENDDO
C$BWD-OF II-LOOP 
      DO j=1,nvor
        DO i=nvor,1,-1
          wc_gam_diff(1, i, j) = wc_gam_diff(1, i, j) + enc(1, i)*
     +      aicn_diff(i, j)
          enc_diff(1, i) = enc_diff(1, i) + wc_gam(1, i, j)*aicn_diff(i
     +      , j)
          wc_gam_diff(2, i, j) = wc_gam_diff(2, i, j) + enc(2, i)*
     +      aicn_diff(i, j)
          enc_diff(2, i) = enc_diff(2, i) + wc_gam(2, i, j)*aicn_diff(i
     +      , j)
          wc_gam_diff(3, i, j) = wc_gam_diff(3, i, j) + enc(3, i)*
     +      aicn_diff(i, j)
          enc_diff(3, i) = enc_diff(3, i) + wc_gam(3, i, j)*aicn_diff(i
     +      , j)
          aicn_diff(i, j) = 0.D0
        ENDDO
      ENDDO
      CALL VVOR_B(betm, iysym, ysym, ysym_diff, izsym, zsym, zsym_diff, 
     +            vrcore, nvor, rv1, rv1_diff, rv2, rv2_diff, nsurfv, 
     +            chordv, chordv_diff, nvor, rc, rc_diff, nsurfv, 
     +            .false., wc_gam, wc_gam_diff, nvmax)
      END

C  Differentiation of velsum in reverse (adjoint) mode (with options i4 dr8 r8):
C   gradient     of useful results: vinf wrot gam wv
C   with respect to varying inputs: vinf wrot gam wv wv_gam
C   RW status of diff variables: vinf:incr wrot:incr gam:incr wv:in-out
C                wv_gam:out
C GAMSUM
C
C
      SUBROUTINE VELSUM_B()
      INCLUDE 'AVL.INC'
      INCLUDE 'AVL_ad_seeds.inc'
      INTEGER i
      INTEGER k
      INTEGER j
      INTEGER n
      INTEGER ii3
      INTEGER ii2
      INTEGER ii1
      DO ii1=1,nvmax
        DO ii2=1,nvmax
          DO ii3=1,3
            wv_gam_diff(ii3, ii2, ii1) = 0.D0
          ENDDO
        ENDDO
      ENDDO
      DO j=nvor,1,-1
        DO i=nvor,1,-1
          DO k=3,1,-1
            wv_gam_diff(k, i, j) = wv_gam_diff(k, i, j) + gam(j)*wv_diff
     +        (k, i)
            gam_diff(j) = gam_diff(j) + wv_gam(k, i, j)*wv_diff(k, i)
          ENDDO
        ENDDO
      ENDDO
      DO i=nvor,1,-1
        DO k=3,1,-1
          vinf_diff(1) = vinf_diff(1) + wvsrd_u(k, i, 1)*wv_diff(k, i)
          vinf_diff(2) = vinf_diff(2) + wvsrd_u(k, i, 2)*wv_diff(k, i)
          vinf_diff(3) = vinf_diff(3) + wvsrd_u(k, i, 3)*wv_diff(k, i)
          wrot_diff(1) = wrot_diff(1) + wvsrd_u(k, i, 4)*wv_diff(k, i)
          wrot_diff(2) = wrot_diff(2) + wvsrd_u(k, i, 5)*wv_diff(k, i)
          wrot_diff(3) = wrot_diff(3) + wvsrd_u(k, i, 6)*wv_diff(k, i)
          wv_diff(k, i) = 0.D0
        ENDDO
      ENDDO
      END

C  Differentiation of set_par_and_cons in reverse (adjoint) mode (with options i4 dr8 r8):
C   gradient     of useful results: alfa beta wrot xyzref
C   with respect to varying inputs: alfa beta conval xyzref
C WSENS
      SUBROUTINE SET_PAR_AND_CONS_B(niter, ir)
      INCLUDE 'AVL.INC'
      INCLUDE 'AVL_ad_seeds.inc'
      INTEGER niter, ir
      INTEGER n
      INTEGER iv
      INTEGER ic
C
C
C
      INTEGER ii2
      INTEGER ii1
      INTEGER branch
      IF (niter .GT. 0) THEN
C----- might as well directly set operating variables if they are known
        IF (icon(ivalfa, ir) .EQ. icalfa) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (icon(ivbeta, ir) .EQ. icbeta) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (icon(ivrotx, ir) .EQ. icrotx) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
        IF (icon(ivroty, ir) .EQ. icroty) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
C$AD-II-loop
        IF (icon(ivrotz, ir) .EQ. icrotz) THEN
          DO ii1=1,nrmax
            DO ii2=1,icmax
              conval_diff(ii2, ii1) = 0.D0
            ENDDO
          ENDDO
          conval_diff(icrotz, ir) = conval_diff(icrotz, ir) + 2.*
     +      wrot_diff(3)/bref
          wrot_diff(3) = 0.D0
        ELSE
          DO ii1=1,nrmax
            DO ii2=1,icmax
              conval_diff(ii2, ii1) = 0.D0
            ENDDO
          ENDDO
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          conval_diff(icroty, ir) = conval_diff(icroty, ir) + 2.*
     +      wrot_diff(2)/cref
          wrot_diff(2) = 0.D0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) conval_diff(icrotx, ir) = conval_diff(icrotx
     +      , ir) + 2.*wrot_diff(1)/bref
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          conval_diff(icbeta, ir) = conval_diff(icbeta, ir) + dtr*
     +      beta_diff
          beta_diff = 0.D0
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          conval_diff(icalfa, ir) = conval_diff(icalfa, ir) + dtr*
     +      alfa_diff
          alfa_diff = 0.D0
        END IF
      ELSE
        DO ii1=1,nrmax
          DO ii2=1,icmax
            conval_diff(ii2, ii1) = 0.D0
          ENDDO
        ENDDO
      END IF
      xyzref_diff(3) = 0.D0
      xyzref_diff(2) = 0.D0
      xyzref_diff(1) = 0.D0
      END

C  Differentiation of set_vel_rhs in reverse (adjoint) mode (with options i4 dr8 r8):
C   gradient     of useful results: vinf xyzref rc rhs
C   with respect to varying inputs: vinf wrot xyzref rc enc
      SUBROUTINE SET_VEL_RHS_B()
C
      INCLUDE 'AVL.INC'
      INCLUDE 'AVL_ad_seeds.inc'
      REAL rrot(3), vunit(3), vunit_w_term(3), wunit(3)
      REAL rrot_diff(3), vunit_diff(3), vunit_w_term_diff(3), wunit_diff
     +     (3)
      INTEGER i
      REAL DOT
      REAL result1
      REAL result1_diff
      INTEGER branch
      DO i=1,nvor
        IF (lvnc(i)) THEN
          CALL PUSHREAL8(vunit(1))
          vunit(1) = 0.
          CALL PUSHREAL8(vunit(2))
          vunit(2) = 0.
          CALL PUSHREAL8(vunit(3))
          vunit(3) = 0.
          CALL PUSHREAL8(wunit(1))
          wunit(1) = 0.
          CALL PUSHREAL8(wunit(2))
          wunit(2) = 0.
          CALL PUSHREAL8(wunit(3))
          wunit(3) = 0.
          IF (lvalbe(i)) THEN
            CALL PUSHREAL8(vunit(1))
            vunit(1) = vinf(1)
            CALL PUSHREAL8(vunit(2))
            vunit(2) = vinf(2)
            CALL PUSHREAL8(vunit(3))
            vunit(3) = vinf(3)
            CALL PUSHREAL8(wunit(1))
            wunit(1) = wrot(1)
            CALL PUSHREAL8(wunit(2))
            wunit(2) = wrot(2)
            CALL PUSHREAL8(wunit(3))
            wunit(3) = wrot(3)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8(rrot(1))
          rrot(1) = rc(1, i) - xyzref(1)
          CALL PUSHREAL8(rrot(2))
          rrot(2) = rc(2, i) - xyzref(2)
          CALL PUSHREAL8(rrot(3))
          rrot(3) = rc(3, i) - xyzref(3)
          CALL CROSS(rrot, wunit, vunit_w_term)
          CALL PUSHREAL8ARRAY(vunit, 3)
          vunit = vunit + vunit_w_term
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      ENDDO
      wrot_diff = 0.D0
      enc_diff = 0.D0
      vunit_diff = 0.D0
      wunit_diff = 0.D0
      vunit_w_term_diff = 0.D0
      rrot_diff = 0.D0
      DO i=nvor,1,-1
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          rhs_diff(i) = 0.D0
        ELSE
          result1_diff = -rhs_diff(i)
          rhs_diff(i) = 0.D0
          CALL DOT_B(enc(1, i), enc_diff(1, i), vunit, vunit_diff, 
     +               result1_diff)
          CALL POPREAL8ARRAY(vunit, 3)
          vunit_w_term_diff = vunit_w_term_diff + vunit_diff
          CALL CROSS_B(rrot, rrot_diff, wunit, wunit_diff, vunit_w_term
     +                 , vunit_w_term_diff)
          CALL POPREAL8(rrot(3))
          rc_diff(3, i) = rc_diff(3, i) + rrot_diff(3)
          xyzref_diff(3) = xyzref_diff(3) - rrot_diff(3)
          rrot_diff(3) = 0.D0
          CALL POPREAL8(rrot(2))
          rc_diff(2, i) = rc_diff(2, i) + rrot_diff(2)
          xyzref_diff(2) = xyzref_diff(2) - rrot_diff(2)
          rrot_diff(2) = 0.D0
          CALL POPREAL8(rrot(1))
          rc_diff(1, i) = rc_diff(1, i) + rrot_diff(1)
          xyzref_diff(1) = xyzref_diff(1) - rrot_diff(1)
          rrot_diff(1) = 0.D0
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(wunit(3))
            wrot_diff(3) = wrot_diff(3) + wunit_diff(3)
            wunit_diff(3) = 0.D0
            CALL POPREAL8(wunit(2))
            wrot_diff(2) = wrot_diff(2) + wunit_diff(2)
            wunit_diff(2) = 0.D0
            CALL POPREAL8(wunit(1))
            wrot_diff(1) = wrot_diff(1) + wunit_diff(1)
            wunit_diff(1) = 0.D0
            CALL POPREAL8(vunit(3))
            vinf_diff(3) = vinf_diff(3) + vunit_diff(3)
            vunit_diff(3) = 0.D0
            CALL POPREAL8(vunit(2))
            vinf_diff(2) = vinf_diff(2) + vunit_diff(2)
            vunit_diff(2) = 0.D0
            CALL POPREAL8(vunit(1))
            vinf_diff(1) = vinf_diff(1) + vunit_diff(1)
            vunit_diff(1) = 0.D0
          END IF
          CALL POPREAL8(wunit(3))
          wunit_diff(3) = 0.D0
          CALL POPREAL8(wunit(2))
          wunit_diff(2) = 0.D0
          CALL POPREAL8(wunit(1))
          wunit_diff(1) = 0.D0
          CALL POPREAL8(vunit(3))
          vunit_diff(3) = 0.D0
          CALL POPREAL8(vunit(2))
          vunit_diff(2) = 0.D0
          CALL POPREAL8(vunit(1))
          vunit_diff(1) = 0.D0
        END IF
      ENDDO
      END

C  Differentiation of mat_prod in reverse (adjoint) mode (with options i4 dr8 r8):
C   gradient     of useful results: vec out_vec
C   with respect to varying inputs: vec mat
Cset_vel_rhs
      SUBROUTINE MAT_PROD_B(mat, mat_diff, vec, vec_diff, n, out_vec, 
     +                      out_vec_diff)
      INCLUDE 'AVL.INC'
      INCLUDE 'AVL_ad_seeds.inc'
      REAL mat(nvmax, nvmax), vec(nvmax), out_vec(nvmax)
      REAL mat_diff(nvmax, nvmax), vec_diff(nvmax), out_vec_diff(nvmax)
      INTEGER j
      INTEGER i
C$AD II-LOOP
      INTEGER n
      mat_diff = 0.D0
C$BWD-OF II-LOOP 
      DO j=1,n
        DO i=n,1,-1
          mat_diff(i, j) = mat_diff(i, j) + vec(j)*out_vec_diff(i)
          vec_diff(j) = vec_diff(j) + mat(i, j)*out_vec_diff(i)
        ENDDO
      ENDDO
      END
Cmat_prod
