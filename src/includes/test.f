      CHARACTER*32 UNCHL,UNCHM,UNCHT,UNCHF,UNCHS,UNCHV,UNCHA,UNCHI,UNCHD
      COMMON /UN_R/UNITL,UNITM,UNITT,UNITF,UNITS,UNITV,UNITA,UNITI,UNITD
      COMMON /UN_C/UNCHL,UNCHM,UNCHT,UNCHF,UNCHS,UNCHV,UNCHA,UNCHI,UNCHD
      COMMON /UN_I/NUL  ,NUM  ,NUT  ,NUF  ,NUS  ,NUV  ,NUA  ,NUI  ,NUD
C
      REAL MACH, MACH0
C
      CHARACTER*80  FRNDEF, FPRDEF, FEVDEF
      CHARACTER*80 TITLE
      CHARACTER*40 STITLE, BTITLE, RTITLE
      CHARACTER*16 DNAME, GNAME
      CHARACTER*12 VARNAM, CONNAM
      CHARACTER*12 VARKEY
      CHARACTER*3  CONKEY
      CHARACTER*10 PARNAM
      CHARACTER*32 PARUNCH
      COMMON /CASE_C/
     & FRNDEF,         ! default run case save file
     & FPRDEF,         ! default dimensional parameter file
     & FEVDEF,         ! default eigenvalue save file
     & TITLE,          ! configuration title
     & STITLE(NFMAX),  ! surface title
     & BTITLE(NBMAX),  ! body title
     & RTITLE(NRMAX),  ! run case title
     & DNAME(NDMAX),   ! control variable name
     & GNAME(NGMAX),   ! design  variable name
     & VARNAM(IVMAX),  ! variable   name
     & CONNAM(ICMAX),  ! constraint name
     & VARKEY(IVMAX),  ! variable   selection key
     & CONKEY(ICMAX),  ! constraint selection key
     & PARNAM(IPMAX),  ! run case parameter name
     & PARUNCH(IPMAX)  ! run case parameter unit name
C
      COMMON /SOLV_R/
     & AMACH,                 ! Mach number at which AIC matrices were computed
     & AICN(NVMAX,NVMAX),     ! normalwash AIC matrix (and VL system matrix)
     & AICN_LU(NVMAX,NVMAX),     ! LU facotrization of AICN
     & WC_GAM(3,NVMAX,NVMAX), ! c.p. velocity/Gamma influence matrix
     & WV_GAM(3,NVMAX,NVMAX), ! h.v. velocity/Gamma influence matrix
     & WC(3,NVMAX),           ! total induced velocity at c.p.
     & WC_U(3,NVMAX,NUMAX),
     & WC_D(3,NVMAX,NDMAX),
     & WC_G(3,NVMAX,NGMAX),
     & WV(3,NVMAX),           ! total induced velocity at h.v.
     & WV_U(3,NVMAX,NUMAX),
     & WV_D(3,NVMAX,NDMAX),
     & WV_G(3,NVMAX,NGMAX)