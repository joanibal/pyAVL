      subroutine write_tecplot(file_name)
C
      INCLUDE 'AVL.INC'
C
      LOGICAL LRANGEALL
      INTEGER LU, nChords, nStrips, nPts
      CHARACTER*(*) file_name
C
C...  check that all surfaces use full range of airfoil definitions
      LRANGEALL = .TRUE.
      DO ISURF = 1, NSURF
        LRANGEALL = LRANGEALL .AND. LRANGE(ISURF)
      ENDDO
      IF (.NOT.LRANGEALL) THEN
      WRITE(*,*) 'ERROR in write_tecplot: implemented only for ',
     &    'surfaces defined using full range of input airfoils'
        WRITE(*,*) '  returning without writing OML output'
        RETURN
      ENDIF
      
c      write out the header
c      we have to use "block packing since the data is all cell cenered"

c      Open the file for writing. Use a unit number (out_unit) that is not in use.
      LU = 13
      OPEN(UNIT=LU, FILE=file_name, STATUS='REPLACE', ACTION='WRITE')
      
      WRITE(LU, '(A,A)') 'TITLE = ', trim(TITLE)
      WRITE(LU, '(A)') 'VARIABLES = "X", "Y", "Z"'
      WRITE(LU, '(A,f5.2)') 'DATASETAUXDATA AVL=',VERSION
      
      WRITE(LU, '(A,A,A)') 'ZONE T=',
     & trim(STITLE(ISURF)),' DATAPACKING=BLOCK'
      
      DO ISURF = 1, NSURF
        nChords = NK(ISURF)
        nStrips = NJ(ISURF)
        
        nPts = (nStrips+1)*(nChords*2+1)
        nPts = (nStrips+1)*(nChords+1)
        DO idim = 1,3
            DO iPt = 1,nPts
                WRITE(LU, '(F15.8)') XYZLO(idim, iPt, ISURF)
            ENDDO 
        ENDDO
        
      ENDDO  
      CLOSE(LU)
      
      end subroutine write_tecplot
      
      
