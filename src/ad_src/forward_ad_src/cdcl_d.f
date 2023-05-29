C        Generated by TAPENADE     (INRIA, Ecuador team)
C  Tapenade 3.16 (develop) - 15 Jan 2021 14:26
C
C  Differentiation of cdcl in forward (tangent) mode (with options i4 dr8 r8):
C   variations   of useful results: cd_cl cd
C   with respect to varying inputs: cl
C***********************************************************************
C    Module:  cdcl.f
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
C
      SUBROUTINE CDCL_D(j, cl, cl_diff, cd, cd_diff, cd_cl, cd_cl_diff)
      INCLUDE 'AVL.INC'
      INCLUDE 'AVL_ad_seeds.inc'
      REAL clmin
      REAL cdmin
      REAL cl0
      REAL cd0
      REAL clmax
      REAL cdmax
      REAL cdx1
      REAL cdx2
      REAL clfac
      REAL temp
      INTEGER j
      REAL clinc
      REAL cd_cl
      REAL cd_cl_diff
      REAL cd
      REAL cd_diff
      REAL cl
      REAL cl_diff
      REAL cdinc
C
C--- CLINC and CDINC control the rate of increase of drag in the stall regions
C    CLINC=0.2 forces drag to increase by CDINC over deltaCL=0.2 after stall
C    The CDINC term is the CD increment added by deltaCL=CLINC after stall
C
      DATA clinc, cdinc /0.2, 0.0500/
C
      cd = 0.0
      cd_cl = 0.0
      IF (j .LT. 1 .OR. j .GT. nstrip) THEN
        WRITE(*, *) 'Error in CDCL - strip index out of bounds'
        cd_cl_diff = 0.D0
        cd_diff = 0.D0
        RETURN
      ELSE
C
C--- Unpack the CD,CL parameters for this strip
        clmin = clcd(1, j)
        cdmin = clcd(2, j)
        cl0 = clcd(3, j)
        cd0 = clcd(4, j)
        clmax = clcd(5, j)
        cdmax = clcd(6, j)
        IF (clmax .LE. cl0 .OR. cl0 .LE. clmin) THEN
          WRITE(*, *) '*** CDCL input CL data out of order'
          cd_cl_diff = 0.D0
          cd_diff = 0.D0
          RETURN
        ELSE
C
C--- Some matching parameters that make the slopes smooth at the stall joins
          cdx1 = 2.0*(cdmin-cd0)*(clmin-cl0)/(clmin-cl0)**2
          cdx2 = 2.0*(cdmax-cd0)*(clmax-cl0)/(clmax-cl0)**2
          clfac = 1.0/clinc
C
C--- Four formulas are used for CD, depending on the CL  
C
          IF (cl .LT. clmin) THEN
C--- Negative stall drag model (slope matches lower side, quadratic drag rise)
            temp = cdinc*(clfac*clfac)
            cd_diff = (temp*2*(cl-clmin)-cdx1/(clmin-cl0))*cl_diff
            cd = cdmin + temp*((cl-clmin)*(cl-clmin)) + cdx1*(1.0-(cl-
     +        cl0)/(clmin-cl0))
            cd_cl_diff = cdinc*2.0*clfac*cl_diff
            cd_cl = cdinc*clfac*2.0*(cl-clmin)
C
C--- Quadratic matching negative stall and minimum drag point with zero slope
          ELSE IF (cl .LT. cl0) THEN
            cd_diff = (cdmin-cd0)*2*(cl-cl0)*cl_diff/(clmin-cl0)**2
            cd = cd0 + (cdmin-cd0)*(cl-cl0)**2/(clmin-cl0)**2
            cd_cl_diff = (cdmin-cd0)*2.0*cl_diff/(clmin-cl0)**2
            cd_cl = (cdmin-cd0)*2.0*(cl-cl0)/(clmin-cl0)**2
C
C--- Quadratic matching positive stall and minimum drag point with zero slope
          ELSE IF (cl .LT. clmax) THEN
            cd_diff = (cdmax-cd0)*2*(cl-cl0)*cl_diff/(clmax-cl0)**2
            cd = cd0 + (cdmax-cd0)*(cl-cl0)**2/(clmax-cl0)**2
            cd_cl_diff = (cdmax-cd0)*2.0*cl_diff/(clmax-cl0)**2
            cd_cl = (cdmax-cd0)*2.0*(cl-cl0)/(clmax-cl0)**2
C
          ELSE
C--- Positive stall drag model (slope matches upper side, quadratic drag rise)
            temp = cdinc*(clfac*clfac)
            cd_diff = (temp*2*(cl-clmax)+cdx2/(clmax-cl0))*cl_diff
            cd = cdmax + temp*((cl-clmax)*(cl-clmax)) - cdx2*(1.0-(cl-
     +        cl0)/(clmax-cl0))
            cd_cl_diff = cdinc*2.0*clfac*cl_diff
            cd_cl = cdinc*clfac*2.0*(cl-clmax)
          END IF
C
          RETURN
        END IF
      END IF
      END
C
C
C
C
C
C
C

