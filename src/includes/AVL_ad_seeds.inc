C
      ! REAL YSYM_diff, ZSYM_diff    ! y- and z-locations of symmetry planes
      ! REAL ALFA_diff, BETA_diff        ! alpha, beta
      ! REAL VINF_diff             ! freestream velocies in body axes
      ! REAL VINF_A_diff           ! d_diff(Vinf)/d_diff(alpha)
      ! REAL VINF_B_diff           ! d_diff(Vinf)/d_diff(beta)
      ! REAL WROT_diff             ! rotation rates in body axes
      ! REAL DTIMED_diff              ! time step for eigenmode movie integration
      ! REAL PARVAL_diff  ! parameter values for run cases
      ! REAL CONVAL_diff  ! imposed constraint value
      ! REAL DELCON_diff        ! imposed control variable value
      ! REAL DELDES_diff        ! imposed design  variable value
      ! REAL SREF_diff, CREF_diff,  BREF_diff   ! Area, Chord, Span  reference values
      ! REAL XYZREF_diff            ! X,Y,Z location for moments
      ! REAL XYZREF0_diff           ! X,Y,Z location for moments _diff(default)
      ! REAL MACH_diff                 ! freestream Mach number
      ! REAL MACH0_diff                ! freestream Mach number _diff(default)
      ! REAL CDREF_diff                ! baseline profile CD
      ! REAL CDREF0_diff               ! baseline profile CD _diff(default)
      ! REAL VRCORE_diff               ! vortex core radius / chord
      ! REAL SRCORE_diff               ! source core radius / body radius
      ! REAL CLFF_diff, CYFF_diff, CDFF_diff     ! Trefftz-plane  CL,CY,CDi
      ! REAL CLFF_U_diff,CYFF_U_diff,CDFF_U_diff ! deriv wrt Vinf,Wrot
      ! REAL CLFF_D_diff,CYFF_D_diff,CDFF_D_diff ! deriv wrt control
      ! REAL CLFF_G_diff,CYFF_G_diff,CDFF_G_diff ! deriv wrt design
      ! REAL SPANEF_diff         ! span efficiency
      ! REAL SPANEF_A_diff       ! d_diff(SPANEF)/d_diff(alpha)
      ! REAL SPANEF_U_diff   ! d_diff(SPANEF)/d_diff(beta)
      ! REAL SPANEF_D_diff   ! d_diff(SPANEF)/d_diff(control)
      ! REAL SPANEF_G_diff   ! d_diff(SPANEF)/d_diff(design)
      ! REAL CDTOT_diff, CLTOT_diff        ! total CD,CL
      ! REAL CXTOT_diff, CYTOT_diff, CZTOT_diff ! total Cx,Cy,Cz
      ! REAL CRTOT_diff, CMTOT_diff, CNTOT_diff ! total Cl,Cm,Cn
      ! REAL CRSAX_diff, CMSAX_diff, CNSAX_diff ! total Cl,Cm,Cn in stability axis
      ! REAL CDVTOT_diff              ! total viscous CD
      ! REAL CDTOT_A_diff, CLTOT_A_diff          ! sensitivities wrt alpha
      ! REAL CDTOT_U_diff,CLTOT_U_diff              !sens wrt U,W
      ! REAL CXTOT_U_diff,CYTOT_U_diff,CZTOT_U_diff!sens wrt U,W
      ! REAL CRTOT_U_diff,CMTOT_U_diff,CNTOT_U_diff!sens wrt U,W
      ! REAL CDTOT_D_diff,CLTOT_D_diff              !sens wrt control
      ! REAL CXTOT_D_diff,CYTOT_D_diff,CZTOT_D_diff!sens wrt control
      ! REAL CRTOT_D_diff,CMTOT_D_diff,CNTOT_D_diff!sens wrt control
      ! REAL CDTOT_G_diff,CLTOT_G_diff              !sens wrt design
      ! REAL CXTOT_G_diff,CYTOT_G_diff,CZTOT_G_diff!sens wrt design
      ! REAL CRTOT_G_diff,CMTOT_G_diff,CNTOT_G_diff!sens wrt design
      ! REAL CHINGE_diff
      ! REAL CHINGE_U_diff
      ! REAL CHINGE_D_diff
      ! REAL CHINGE_G_diff
      ! REAL DCL_A0_diff, DCM_A0_diff ! additional default CL_a, CM_a
      ! REAL DCL_U0_diff, DCM_U0_diff  ! additional default CL_u, CM_u
      ! REAL XNP_diff                ! x neutral point
      
      COMMON /CASE_R_d/
     & YSYM_diff,   ZSYM_diff,    ! y- and z-locations of symmetry planes
     & ALFA_diff,   BETA_diff,        ! alpha, beta
     & VINF_diff(3),             ! freestream velocies in body axes
     & VINF_A_diff(3),           ! d_diff(Vinf)/d_diff(alpha)
     & VINF_B_diff(3),           ! d_diff(Vinf)/d_diff(beta)
     & WROT_diff(3),             ! rotation rates in body axes
     & DTIMED_diff,              ! time step for eigenmode movie integration
     & PARVAL_diff(IPMAX,NRMAX),  ! parameter values for run cases
     & CONVAL_diff(ICMAX,NRMAX),  ! imposed constraint value
     & DELCON_diff(NDMAX),        ! imposed control variable value
     & DELDES_diff(NGMAX),        ! imposed design  variable value
     & SREF_diff,  CREF_diff,  BREF_diff,   ! Area, Chord, Span  reference values
     & XYZREF_diff(3),            ! X,Y,Z location for moments
     & XYZREF0_diff(3),           ! X,Y,Z location for moments _diff(default)
     & MACH_diff,                 ! freestream Mach number
     & MACH0_diff,                ! freestream Mach number _diff(default)
     & CDREF_diff,                ! baseline profile CD
     & CDREF0_diff,               ! baseline profile CD _diff(default)
     & VRCORE_diff,               ! vortex core radius / chord
     & SRCORE_diff,               ! source core radius / body radius
     & CLFF_diff, CYFF_diff, CDFF_diff,     ! Trefftz-plane  CL,CY,CDi
     & CLFF_U_diff(NUMAX),CYFF_U_diff(NUMAX),CDFF_U_diff(NUMAX), ! deriv wrt Vinf,Wrot
     & CLFF_D_diff(NDMAX),CYFF_D_diff(NDMAX),CDFF_D_diff(NDMAX), ! deriv wrt control
     & CLFF_G_diff(NGMAX),CYFF_G_diff(NGMAX),CDFF_G_diff(NGMAX), ! deriv wrt design
     & SPANEF_diff,           ! span efficiency
     & SPANEF_A_diff,         ! d_diff(SPANEF)/d_diff(alpha)
     & SPANEF_U_diff(NUMAX),  ! d_diff(SPANEF)/d_diff(beta)
     & SPANEF_D_diff(NDMAX),  ! d_diff(SPANEF)/d_diff(control)
     & SPANEF_G_diff(NGMAX),  ! d_diff(SPANEF)/d_diff(design)
     & CDTOT_diff, CLTOT_diff,          ! total CD,CL
     & CXTOT_diff, CYTOT_diff, CZTOT_diff,   ! total Cx,Cy,Cz
     & CRTOT_diff, CMTOT_diff, CNTOT_diff,   ! total Cl,Cm,Cn
     & CRSAX_diff, CMSAX_diff, CNSAX_diff,   ! total Cl,Cm,Cn in stability axis
     & CDVTOT_diff,                ! total viscous CD
     & CDTOT_A_diff, CLTOT_A_diff,            ! sensitivities wrt alpha
     & CDTOT_U_diff(NUMAX),CLTOT_U_diff(NUMAX),                !sens wrt U,W
     & CXTOT_U_diff(NUMAX),CYTOT_U_diff(NUMAX),CZTOT_U_diff(NUMAX), !sens wrt U,W
     & CRTOT_U_diff(NUMAX),CMTOT_U_diff(NUMAX),CNTOT_U_diff(NUMAX), !sens wrt U,W
     & CDTOT_D_diff(NDMAX),CLTOT_D_diff(NDMAX),                !sens wrt control
     & CXTOT_D_diff(NDMAX),CYTOT_D_diff(NDMAX),CZTOT_D_diff(NDMAX), !sens wrt control
     & CRTOT_D_diff(NDMAX),CMTOT_D_diff(NDMAX),CNTOT_D_diff(NDMAX), !sens wrt control
     & CDTOT_G_diff(NGMAX),CLTOT_G_diff(NGMAX),                !sens wrt design
     & CXTOT_G_diff(NGMAX),CYTOT_G_diff(NGMAX),CZTOT_G_diff(NGMAX), !sens wrt design
     & CRTOT_G_diff(NGMAX),CMTOT_G_diff(NGMAX),CNTOT_G_diff(NGMAX), !sens wrt design
     & CHINGE_diff(NDMAX),
     & CHINGE_U_diff(NDMAX,NUMAX),
     & CHINGE_D_diff(NDMAX,NDMAX),
     & CHINGE_G_diff(NDMAX,NGMAX),
     & DCL_A0_diff, DCM_A0_diff,   ! additional default CL_a, CM_a
     & DCL_U0_diff, DCM_U0_diff,    ! additional default CL_u, CM_u
     & XNP_diff                ! x neutral point

      COMMON /SURF_R_d/
     & CDSURF_diff(NFMAX),CLSURF_diff(NFMAX),                ! surface CD,CL
     & CXSURF_diff(NFMAX),CYSURF_diff(NFMAX),CZSURF_diff(NFMAX),  ! surface Cx,Cy,Cz
     & CRSURF_diff(NFMAX),CNSURF_diff(NFMAX),CMSURF_diff(NFMAX),  ! surface Cl,Cm,Cn
     & CDVSURF_diff(NFMAX),                             ! surface viscous CD
     & CDS_A_diff(NFMAX), CLS_A_diff(NFMAX),
     & CDS_U_diff(NFMAX,NUMAX),
     & CLS_U_diff(NFMAX,NUMAX),
     & CXS_U_diff(NFMAX,NUMAX),
     & CYS_U_diff(NFMAX,NUMAX),
     & CZS_U_diff(NFMAX,NUMAX),
     & CRS_U_diff(NFMAX,NUMAX),
     & CNS_U_diff(NFMAX,NUMAX),
     & CMS_U_diff(NFMAX,NUMAX),
     & CDS_D_diff(NFMAX,NDMAX),
     & CLS_D_diff(NFMAX,NDMAX),
     & CXS_D_diff(NFMAX,NDMAX),
     & CYS_D_diff(NFMAX,NDMAX),
     & CZS_D_diff(NFMAX,NDMAX),
     & CRS_D_diff(NFMAX,NDMAX),
     & CNS_D_diff(NFMAX,NDMAX),
     & CMS_D_diff(NFMAX,NDMAX),
     & CDS_G_diff(NFMAX,NGMAX),
     & CLS_G_diff(NFMAX,NGMAX),
     & CXS_G_diff(NFMAX,NGMAX),
     & CYS_G_diff(NFMAX,NGMAX),
     & CZS_G_diff(NFMAX,NGMAX),
     & CRS_G_diff(NFMAX,NGMAX),
     & CNS_G_diff(NFMAX,NGMAX),
     & CMS_G_diff(NFMAX,NGMAX),
     & CF_SRF_diff(3,NFMAX), CM_SRF_diff(3,NFMAX),
     & CL_SRF_diff(NFMAX), CD_SRF_diff(NFMAX), CMLE_SRF_diff(NFMAX),
     & SSURF_diff(NFMAX), CAVESURF_diff(NFMAX),
     & VBODY_diff(NBMAX)
c
c
      COMMON /SURF_GEOM_R_d/
     & XYZSCAL_diff(3, NFMAX), ! scaling factors for XYZ coordinates
     & XYZTRAN_diff(3, NFMAX), ! translation factors for XYZ coordinates
     & YDUPL_diff(NFMAX), ! y position of duplicate plane
     & ADDINC_diff(NFMAX), ! additional incidence angle
     & CSPACE_diff(NFMAX), ! chordwise spacing
     & SSPACE_diff(NFMAX), ! spanwise spacing
     & SSPACES_diff(NSMAX, NFMAX), ! spanwise spacing vector
     & XYZLES_diff(3, NSMAX, NFMAX), ! leading edge cordinate vector
     & CHORDS_diff(NSMAX, NFMAX), ! chord length vectorm
     & AINCS_diff(NSMAX, NFMAX), ! incidence angle vector
     & XASEC_diff(IBX, NSMAX, NFMAX),  ! the x coordinate aifoil section
     & SASEC_diff(IBX, NSMAX, NFMAX),  ! slope of camber line at xasec
     & TASEC_diff(IBX, NSMAX, NFMAX),   ! thickness at xasec
     & CLCDSEC_diff(6, NSMAX, NFMAX),
     & CLAF_diff(NSMAX, NFMAX),
     & XHINGED_diff(ICONX, NSMAX, NFMAX), ! hinge location
     & VHINGED_diff(3, ICONX, NSMAX, NFMAX), ! hinge vector
     & GAIND_diff(ICONX, NSMAX, NFMAX), ! control surface gain
     & REFLD_diff(ICONX, NSMAX, NFMAX), ! control surface reflection
     & GAING_diff(ICONX, NSMAX, NFMax)  ! desgin variable gain
c
c
      COMMON /STRP_R_d/
     & RLE_diff(3,NSMAX),  CHORD_diff(NSMAX),    ! strip c.p. line LE point, chord
     & RLE1_diff(3,NSMAX), CHORD1_diff(NSMAX),   ! strip left  end LE point, chord
     & RLE2_diff(3,NSMAX), CHORD2_diff(NSMAX),   ! strip right end LE point, chord
     & WSTRIP_diff(NSMAX),                  ! strip y-z width
     & TANLE_diff(NSMAX),  TANTE_diff(NSMAX),    ! strip LE,TE sweep slopes
     & CLCD_diff(NUMAX,NSMAX),              ! strip viscous polar
     & SAXFR_diff,                       ! x/c of spanwise axis for Vperp def
     & ESS_diff(3,NSMAX),                ! spanwise unit vector for Vperp def
     & ENSY_diff(NSMAX), ENSZ_diff(NSMAX),    ! strip normal vector in Trefftz-Pln
     & XSREF_diff(NSMAX),YSREF_diff(NSMAX),ZSREF_diff(NSMAX),
     & AINC_diff(NSMAX),              ! strip's incidence twist angle
     & AINC_G_diff(NSMAX,NGMAX),      ! dAINC/dG
     & CDSTRP_diff(NSMAX),   CLSTRP_diff(NSMAX),                   ! strip forces
     & CXSTRP_diff(NSMAX),   CYSTRP_diff(NSMAX),   CZSTRP_diff(NSMAX),  ! strip forces
     & CRSTRP_diff(NSMAX),   CNSTRP_diff(NSMAX),   CMSTRP_diff(NSMAX),  ! strip moments
     & CDST_A_diff(NSMAX),   CLST_A_diff(NSMAX),                   ! alpha sens.
     & CDST_U_diff(NSMAX,NUMAX),
     & CLST_U_diff(NSMAX,NUMAX),
     & CXST_U_diff(NSMAX,NUMAX),
     & CYST_U_diff(NSMAX,NUMAX),
     & CZST_U_diff(NSMAX,NUMAX),
     & CRST_U_diff(NSMAX,NUMAX),
     & CNST_U_diff(NSMAX,NUMAX),
     & CMST_U_diff(NSMAX,NUMAX),
     & CDST_D_diff(NSMAX,NDMAX),
     & CLST_D_diff(NSMAX,NDMAX),
     & CXST_D_diff(NSMAX,NDMAX),
     & CYST_D_diff(NSMAX,NDMAX),
     & CZST_D_diff(NSMAX,NDMAX),
     & CRST_D_diff(NSMAX,NDMAX),
     & CNST_D_diff(NSMAX,NDMAX),
     & CMST_D_diff(NSMAX,NDMAX),
     & CDST_G_diff(NSMAX,NGMAX),
     & CLST_G_diff(NSMAX,NGMAX),
     & CXST_G_diff(NSMAX,NGMAX),
     & CYST_G_diff(NSMAX,NGMAX),
     & CZST_G_diff(NSMAX,NGMAX),
     & CRST_G_diff(NSMAX,NGMAX),
     & CNST_G_diff(NSMAX,NGMAX),
     & CMST_G_diff(NSMAX,NGMAX),
     & CF_STRP_diff(3,NSMAX),   CM_STRP_diff(3,NSMAX), ! strip forces, norm local
     & CNRMSTRP_diff(NSMAX),    CAXLSTRP_diff(NSMAX),
     & CD_LSTRP_diff(NSMAX),  CL_LSTRP_diff(NSMAX),
     & CDV_LSTRP_diff(NSMAX), CLTSTRP_diff(NSMAX), CLASTRP_diff(NSMAX),
     & CMC4_diff(NSMAX),        CMLE_diff(NSMAX),
     & CNC_diff(NSMAX),         DWWAKE_diff(NSMAX),
     & CNC_U_diff(NSMAX,NUMAX),
     & CNC_D_diff(NSMAX,NDMAX),
     & CNC_G_diff(NSMAX,NGMAX)
C
      COMMON /VRTX_R_d/
     & RV1_diff(3,NVMAX),    ! h.v. vortex left  points
     & RV2_diff(3,NVMAX),    ! h.v. vortex right points
     & RV_diff(3,NVMAX),     ! h.v. vortex center points
     & RC_diff(3,NVMAX),     ! h.v. control points
     & RS_diff(3,NVMAX),     ! h.v. source points
     & RL_diff(3,NLMAX),RADL_diff(NLMAX),  ! source line node points, body radius
     & DXV_diff(NVMAX),               ! chord of element
     & CHORDV_diff(NVMAX),            ! chord of element-containing strip
     & SLOPEV_diff(NVMAX),            ! camber slopes at h.v. bound leg
     & SLOPEC_diff(NVMAX),            ! camber slopes at c.p.
     & DCONTROL_diff(NVMAX,NDMAX),  ! d_diff(normal angle)/dCONTROL
     & VHINGE_diff(3,NSMAX,NDMAX),  ! hinge vector for CONTROL rot. of normal
     & PHINGE_diff(3,NSMAX,NDMAX),  ! point on hingeline for hinge moment calculation
     & VREFL_diff(NSMAX,NDMAX),  !sign applied to hinge vec. of refl. surface
     & ENC_diff(3,NVMAX),           ! c.p. normal vector
     & ENV_diff(3,NVMAX),           ! h.v. normal vector
     & ENC_D_diff(3,NVMAX,NDMAX),   ! sensitivities
     & ENC_G_diff(3,NVMAX,NGMAX),
     & ENV_D_diff(3,NVMAX,NDMAX),
     & ENV_G_diff(3,NVMAX,NGMAX),
     & DCP_diff(NVMAX),             ! delta_diff(Cp) on vortex element
     & DCP_U_diff(NVMAX,NUMAX),
     & DCP_D_diff(NVMAX,NDMAX),
     & DCP_G_diff(NVMAX,NGMAX),
     & GAM_diff(NVMAX),             ! circulation of h.v. vortex
     & GAM_U_0_diff(NVMAX,NUMAX),
     & GAM_U_D_diff(NVMAX,NUMAX,NDMAX),
     & GAM_U_G_diff(NVMAX,NUMAX,NGMAX),
     & GAM_U_diff(NVMAX,NUMAX),
     & GAM_D_diff(NVMAX,NDMAX),
     & GAM_G_diff(NVMAX,NGMAX),
     & SRC_diff(NLMAX),        ! source  strength of source+doublet line elem
     & DBL_diff(3,NLMAX),      ! doublet strength of source+doublet line elem
     & SRC_U_diff(NLMAX,NUMAX),
     & DBL_U_diff(3,NLMAX,NUMAX),
     & WCSRD_U_diff(3,NVMAX,NUMAX),
     & WVSRD_U_diff(3,NVMAX,NUMAX),
     & RHS_diff(NVMAX),             ! RHS _diff(flow normal to surface) for linear system
     & RES_diff(NVMAX)             ! residual for linear system res = Ax-b      
C
      COMMON /BODY_R_d/
     & ELBDY_diff(NBMAX),  ! body length
     & SRFBDY_diff(NBMAX), ! body surface area
     & VOLBDY_diff(NBMAX), ! body volume
     & DCPB_diff(3,NLMAX),                           ! dCp loading on sources
     & CDBDY_diff(NBMAX),CLBDY_diff(NBMAX),               ! body CD,CL
     & CXBDY_diff(NBMAX),CYBDY_diff(NBMAX),CZBDY_diff(NBMAX),  ! body CX,CY,CZ
     & CRBDY_diff(NBMAX),CNBDY_diff(NBMAX),CMBDY_diff(NBMAX)   ! body Cl,Cm,Cn

      COMMON /SOLV_R_d/
     & AMACH_diff,                 ! Mach number at which AIC matrices were computed
     & AICN_diff(NVMAX,NVMAX),     ! normalwash AIC matrix _diff(and VL system matrix)
     & AICN_LU_diff(NVMAX,NVMAX),     ! LU facotrization of AICN
     & WC_GAM_diff(3,NVMAX,NVMAX), ! c.p. velocity/Gamma influence matrix
     & WV_GAM_diff(3,NVMAX,NVMAX), ! h.v. velocity/Gamma influence matrix
     & WC_diff(3,NVMAX),           ! total induced velocity at c.p.
     & WC_U_diff(3,NVMAX,NUMAX),
     & WC_D_diff(3,NVMAX,NDMAX),
     & WC_G_diff(3,NVMAX,NGMAX),
     & WV_diff(3,NVMAX),           ! total induced velocity at h.v.
     & WV_U_diff(3,NVMAX,NUMAX),
     & WV_D_diff(3,NVMAX,NDMAX),
     & WV_G_diff(3,NVMAX,NGMAX)     