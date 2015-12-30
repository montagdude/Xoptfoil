C  This file is part of XOPTFOIL

C  XOPTFOIL is free software: you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation, either version 3 of the License, or
C  (at your option) any later version.

C  XOPTFOIL is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  You should have received a copy of the GNU General Public License
C  along with XOPTFOIL.  If not, see <http://www.gnu.org/licenses/>.

C  Copyright (C) 2014 -- 2015 Daniel Prosser (this modified version of 
C  XFoil code)
C  Original copyright (C) 2000 Mark Drela (original XFoil code)

C===================================================================70
C     *******************************************************
C     *                                                     *
C     *   Factors a full NxN matrix into an LU form.        *
C     *   Subr. BAKSUB can back-substitute it with some RHS.*
C     *   Assumes matrix is non-singular...                 *
C     *    ...if it isn't, a divide by zero will result.    *
C     *                                                     *
C     *   A is the matrix...                                *
C     *     ...replaced with its LU factors.                *
C     *                                                     *
C     *                              Mark Drela  1988       *
C     *******************************************************
C===================================================================70
      SUBROUTINE LUDCMP(NSIZ,N,A,INDX)
C
      DIMENSION A(NSIZ,NSIZ), INDX(NSIZ)
C
      PARAMETER (NVX=500)
      DIMENSION VV(NVX)
C
      IF(N.GT.NVX) STOP 'LUDCMP: Array overflow. Increase NVX.'
C
      DO 12 I=1, N
        AAMAX = 0.
        DO 11 J=1, N
          AAMAX = MAX( ABS(A(I,J)) , AAMAX )
   11   CONTINUE
        VV(I) = 1.0/AAMAX
   12 CONTINUE
C
      DO 19 J=1, N
        DO 14 I=1, J-1
          SUM = A(I,J)
          DO 13 K=1, I-1
            SUM = SUM - A(I,K)*A(K,J)
   13     CONTINUE
          A(I,J) = SUM
   14   CONTINUE
C
        AAMAX = 0.
        DO 16 I=J, N
          SUM = A(I,J)
          DO 15 K=1, J-1
            SUM = SUM - A(I,K)*A(K,J)
   15     CONTINUE
          A(I,J) = SUM
C
          DUM = VV(I)*ABS(SUM)
          IF(DUM.GE.AAMAX) THEN
           IMAX = I
           AAMAX = DUM
          ENDIF
   16   CONTINUE
C
        IF(J.NE.IMAX) THEN
         DO 17 K=1, N
           DUM = A(IMAX,K)
           A(IMAX,K) = A(J,K)
           A(J,K) = DUM
   17    CONTINUE
         VV(IMAX) = VV(J)
        ENDIF
C
        INDX(J) = IMAX
        IF(J.NE.N) THEN
         DUM = 1.0/A(J,J)
         DO 18 I=J+1, N
           A(I,J) = A(I,J)*DUM
   18    CONTINUE
        ENDIF
C
   19 CONTINUE
C
      RETURN
      END ! LUDCMP


C===================================================================70
C===================================================================70
      SUBROUTINE BAKSUB(NSIZ,N,A,INDX,B)
      DIMENSION A(NSIZ,NSIZ), B(NSIZ), INDX(NSIZ)
C
      II = 0
      DO 12 I=1, N
        LL = INDX(I)
        SUM = B(LL)
        B(LL) = B(I)
        IF(II.NE.0) THEN
         DO 11 J=II, I-1
           SUM = SUM - A(I,J)*B(J)
   11    CONTINUE
        ELSE IF(SUM.NE.0.0) THEN
         II = I
        ENDIF
        B(I) = SUM
   12 CONTINUE
C
      DO 14 I=N, 1, -1
        SUM = B(I)
        IF(I.LT.N) THEN
         DO 13 J=I+1, N
           SUM = SUM - A(I,J)*B(J)
   13    CONTINUE
        ENDIF
        B(I) = SUM/A(I,I)
   14 CONTINUE
C
      RETURN
      END ! BAKSUB

C===================================================================70
C
C     Solves KK long, tri-diagonal system |
C                                         |
C             A C          D              |
C             B A C        D              |
C               B A .      .              |
C                 . . C    .              |
C                   B A    D              |
C                                         |
C     The righthand side D is replaced by |
C     the solution.  A, C are destroyed.  |
C
C===================================================================70
      SUBROUTINE TRISOL(A,B,C,D,KK)
      DIMENSION A(KK),B(KK),C(KK),D(KK)
C
      DO 1 K=2, KK
        KM = K-1
        C(KM) = C(KM) / A(KM)
        D(KM) = D(KM) / A(KM)
        A(K) = A(K) - B(K)*C(KM)
        D(K) = D(K) - B(K)*D(KM)
    1 CONTINUE
C
      D(KK) = D(KK)/A(KK)
C
      DO 2 K=KK-1, 1, -1
        D(K) = D(K) - C(K)*D(K+1)
    2 CONTINUE
C
      RETURN
      END ! TRISOL

C===================================================================70
C---------------------------------------------------------------
C     ATAN2 function with branch cut checking.
C
C     Increments position angle of point X,Y from some previous
C     value THOLD due to a change in position, ensuring that the
C     position change does not cross the ATAN2 branch cut
C     (which is in the -x direction).  For example:
C
C       ATANC( -1.0 , -1.0 , 0.75*pi )  returns  1.25*pi , whereas
C       ATAN2( -1.0 , -1.0 )            returns  -.75*pi .
C
C     Typically, ATANC is used to fill an array of angles:
C
C        THETA(1) = ATAN2( Y(1) , X(1) )
C        DO i=2, N
C          THETA(i) = ATANC( Y(i) , X(i) , THETA(i-1) )
C        END DO
C
C     This will prevent the angle array THETA(i) from jumping by 
C     +/- 2 pi when the path X(i),Y(i) crosses the negative x axis.
C
C     Input:
C       X,Y     point position coordinates
C       THOLD   position angle of nearby point
C
C     Output:
C       ATANC   position angle of X,Y
C---------------------------------------------------------------
C===================================================================70
      FUNCTION ATANC(Y,X,THOLD)
      IMPLICIT REAL (A-H,M,O-Z)
      DATA  PI /3.1415926535897932384/
      DATA TPI /6.2831853071795864769/
C
C---- set new position angle, ignoring branch cut in ATAN2 function for now
      THNEW = ATAN2( Y , X )
      DTHET = THNEW - THOLD
C
C---- angle change cannot exceed +/- pi, so get rid of any multiples of 2 pi 
      DTCORR = DTHET - TPI*INT( (DTHET + SIGN(PI,DTHET))/TPI )
C
C---- set correct new angle
      ATANC = THOLD + DTCORR
C
      RETURN
      END ! ATANC

C===================================================================70
C     *******************************************************
C     *                                                     *
C     *   Solves general NxN system in NN unknowns          *
C     *    with arbitrary number (NRHS) of righthand sides. *
C     *   Assumes system is invertible...                   *
C     *    ...if it isn't, a divide by zero will result.    *
C     *                                                     *
C     *   Z is the coefficient matrix...                    *
C     *     ...destroyed during solution process.           *
C     *   R is the righthand side(s)...                     *
C     *     ...replaced by the solution vector(s).          *
C     *                                                     *
C     *                              Mark Drela  1984       *
C     *******************************************************
C===================================================================70
      SUBROUTINE GAUSS(NSIZ,NN,Z,R,NRHS)
C
      DIMENSION Z(NSIZ,NSIZ), R(NSIZ,NRHS)
C
      DO 1 NP=1, NN-1
        NP1 = NP+1
C
C------ find max pivot index NX
        NX = NP
        DO 11 N=NP1, NN
          IF(ABS(Z(N,NP))-ABS(Z(NX,NP))) 11,11,111
  111      NX = N
   11   CONTINUE
C
        PIVOT = 1.0/Z(NX,NP)
C
C------ switch pivots
        Z(NX,NP) = Z(NP,NP)
C
C------ switch rows & normalize pivot row
        DO 12 L=NP1, NN
          TEMP = Z(NX,L)*PIVOT
          Z(NX,L) = Z(NP,L)
          Z(NP,L) = TEMP
   12   CONTINUE
C
        DO 13 L=1, NRHS
          TEMP = R(NX,L)*PIVOT
          R(NX,L) = R(NP,L)
          R(NP,L) = TEMP
   13   CONTINUE
C
C------ forward eliminate everything
        DO 15 K=NP1, NN
          ZTMP = Z(K,NP)
C
C          IF(ZTMP.EQ.0.0) GO TO 15
C
          DO 151 L=NP1, NN
            Z(K,L) = Z(K,L) - ZTMP*Z(NP,L)
  151     CONTINUE
          DO 152 L=1, NRHS
            R(K,L) = R(K,L) - ZTMP*R(NP,L)
  152     CONTINUE
   15   CONTINUE
C
    1 CONTINUE
C
C---- solve for last row
      DO 2 L=1, NRHS
        R(NN,L) = R(NN,L)/Z(NN,NN)
    2 CONTINUE
C
C---- back substitute everything
      DO 3 NP=NN-1, 1, -1
        NP1 = NP+1
        DO 31 L=1, NRHS
          DO 310 K=NP1, NN
            R(NP,L) = R(NP,L) - Z(NP,K)*R(K,L)
  310     CONTINUE
   31   CONTINUE
    3 CONTINUE
C
      RETURN
      END ! GAUSS
