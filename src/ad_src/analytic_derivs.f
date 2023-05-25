      SUBROUTINE solve_adjoint():
        IMPLICIT NONE
        INCLUDE 'AVL_d.inc'
        ! REAL WORK(NVMAX)
        REAL df_dr(NVOR)
        
        CALL BAKSUB(NVMAX,NVOR,AICN_LU,IAPIV,df_dr)
        
        
        
      END ! ADJOINT
        
        
        