C  This file is part of XOPTFOIL
C       Note: not needed, because BIJ is not used again after this

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
C
C     Calculates total and projected TE gap 
C     areas and TE panel strengths.
C
C===================================================================70
      SUBROUTINE TECALC

      use xfoil_inc
C
C---- set TE base vector and TE bisector components
      DXTE = X(1) - X(N)
      DYTE = Y(1) - Y(N)
      DXS = 0.5*(-XP(1) + XP(N))
      DYS = 0.5*(-YP(1) + YP(N))
C
C---- normal and streamwise projected TE gap areas
      ANTE = DXS*DYTE - DYS*DXTE
      ASTE = DXS*DXTE + DYS*DYTE
C
C---- total TE gap area
      DSTE = SQRT(DXTE**2 + DYTE**2)
C
      SHARP = DSTE .LT. 0.0001*CHORD
C
      IF(SHARP) THEN
       SCS = 1.0
       SDS = 0.0
      ELSE
       SCS = ANTE/DSTE
       SDS = ASTE/DSTE
      ENDIF
C
C---- TE panel source and vorticity strengths
      SIGTE = 0.5*(GAM(1) - GAM(N))*SCS
      GAMTE = -.5*(GAM(1) - GAM(N))*SDS
C
      SIGTE_A = 0.5*(GAM_A(1) - GAM_A(N))*SCS
      GAMTE_A = -.5*(GAM_A(1) - GAM_A(N))*SDS
C
      RETURN
      END ! TECALC

C===================================================================70
C
C     Calculates the arc length array S  |
C     for a 2-D array of points (X,Y).   |
C
C===================================================================70
      SUBROUTINE SCALC(X,Y,S,N)
      DIMENSION X(N), Y(N), S(N)
C
      S(1) = 0.
      DO 10 I=2, N
        S(I) = S(I-1) + SQRT((X(I)-X(I-1))**2 + (Y(I)-Y(I-1))**2)
   10 CONTINUE
C
      RETURN
      END ! SCALC

C===================================================================70
C
C     Calculates spline coefficients for X(S).          |
C     Specified 1st derivative and/or usual zero 2nd    |
C     derivative end conditions are used.               |
C     To evaluate the spline at some value of S,        |
C     use SEVAL and/or DEVAL.                           |
C                                                       |
C     S        independent variable array (input)       |
C     X        dependent variable array   (input)       |
C     XS       dX/dS array                (calculated)  |
C     N        number of points           (input)       |
C     XS1,XS2  endpoint derivatives       (input)       |
C              If = 999.0, then usual zero second       |
C              derivative end condition(s) are used     |
C              If = -999.0, then zero third             |
C              derivative end condition(s) are used     |
C                                                       |
C
C===================================================================70
      SUBROUTINE SPLIND(X,XS,S,N,XS1,XS2)
      DIMENSION X(N),XS(N),S(N)
      PARAMETER (NMAX=1000)
      DIMENSION  A(NMAX),B(NMAX),C(NMAX)
      IF(N.GT.NMAX) STOP 'SPLIND: array overflow, increase NMAX'
C     
      DO 1 I=2, N-1
        DSM = S(I) - S(I-1)
        DSP = S(I+1) - S(I)
        B(I) = DSP
        A(I) = 2.0*(DSM+DSP)
        C(I) = DSM
        XS(I) = 3.0*((X(I+1)-X(I))*DSM/DSP + (X(I)-X(I-1))*DSP/DSM)
    1 CONTINUE
C
      IF(XS1.EQ.999.0) THEN
C----- set zero second derivative end condition
       A(1) = 2.0
       C(1) = 1.0
       XS(1) = 3.0*(X(2)-X(1)) / (S(2)-S(1))
      ELSE IF(XS1.EQ.-999.0) THEN
C----- set zero third derivative end condition
       A(1) = 1.0
       C(1) = 1.0
       XS(1) = 2.0*(X(2)-X(1)) / (S(2)-S(1))
      ELSE
C----- set specified first derivative end condition
       A(1) = 1.0
       C(1) = 0.
       XS(1) = XS1
      ENDIF
C
      IF(XS2.EQ.999.0) THEN
       B(N) = 1.0
       A(N) = 2.0
       XS(N) = 3.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
      ELSE IF(XS2.EQ.-999.0) THEN
       B(N) = 1.0
       A(N) = 1.0
       XS(N) = 2.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
      ELSE
       A(N) = 1.0
       B(N) = 0.
       XS(N) = XS2
      ENDIF
C
      IF(N.EQ.2 .AND. XS1.EQ.-999.0 .AND. XS2.EQ.-999.0) THEN
       B(N) = 1.0
       A(N) = 2.0
       XS(N) = 3.0*(X(N)-X(N-1)) / (S(N)-S(N-1))
      ENDIF
C
C---- solve for derivative array XS
      CALL TRISOL(A,B,C,XS,N)
C
      RETURN
      END ! SPLIND

C===================================================================70
C
C     Splines X(S) array just like SPLINE,      |
C     but allows derivative discontinuities     |
C     at segment joints.  Segment joints are    |
C     defined by identical successive S values. |
C
C===================================================================70
      SUBROUTINE SEGSPL(X,XS,S,N)
      DIMENSION X(N), XS(N), S(N)
C
      IF(S(1).EQ.S(2)  ) STOP 'SEGSPL:  First input point duplicated'
      IF(S(N).EQ.S(N-1)) STOP 'SEGSPL:  Last  input point duplicated'
C
      ISEG0 = 1
      DO 10 ISEG=2, N-2
        IF(S(ISEG).EQ.S(ISEG+1)) THEN
         NSEG = ISEG - ISEG0 + 1
         CALL SPLIND(X(ISEG0),XS(ISEG0),S(ISEG0),NSEG,-999.0,-999.0)
         ISEG0 = ISEG+1
        ENDIF
   10 CONTINUE
C
      NSEG = N - ISEG0 + 1
      CALL SPLIND(X(ISEG0),XS(ISEG0),S(ISEG0),NSEG,-999.0,-999.0)
C
      RETURN
      END ! SEGSPL

C===================================================================70
C
C     Calculates dX/dS(SS)                         |
C     XS array must have been calculated by SPLINE |
C
C===================================================================70
      FUNCTION DEVAL(SS,X,XS,S,N)
      DIMENSION X(N), XS(N), S(N)
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      DEVAL = X(I) - X(I-1) + (1.-4.0*T+3.0*T*T)*CX1 + T*(3.0*T-2.)*CX2
      DEVAL = DEVAL/DS
      RETURN
      END ! DEVAL

C===================================================================70
C
C     Calculates d2X/dS2(SS)                       |
C     XS array must have been calculated by SPLINE |
C
C===================================================================70
      FUNCTION D2VAL(SS,X,XS,S,N)
      DIMENSION X(N), XS(N), S(N)

      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      D2VAL = (6.*T-4.)*CX1 + (6.*T-2.0)*CX2
      D2VAL = D2VAL/DS**2
      RETURN
      END ! D2VAL

C===================================================================70
C
C     Locates leading edge spline-parameter value SLE
C
C     The defining condition is
C         
C      (X-XTE,Y-YTE) . (X',Y') = 0     at  S = SLE
C
C     i.e. the surface tangent is normal to the chord
C     line connecting X(SLE),Y(SLE) and the TE point.
C
C===================================================================70
      SUBROUTINE LEFIND(SLE,X,XP,Y,YP,S,N,SILENT_MODE)
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
      LOGICAL SILENT_MODE
C
C---- convergence tolerance
      DSEPS = (S(N)-S(1)) * 1.0E-5
C
C---- set trailing edge point coordinates
      XTE = 0.5*(X(1) + X(N))
      YTE = 0.5*(Y(1) + Y(N))
C
C---- get first guess for SLE
      DO 10 I=3, N-2
        DXTE = X(I) - XTE
        DYTE = Y(I) - YTE
        DX = X(I+1) - X(I)
        DY = Y(I+1) - Y(I)
        DOTP = DXTE*DX + DYTE*DY
        IF(DOTP .LT. 0.0) GO TO 11
   10 CONTINUE
C
   11 SLE = S(I)
C
C---- check for sharp LE case
      IF(S(I) .EQ. S(I-1)) THEN
ccc        WRITE(*,*) 'Sharp LE found at ',I,SLE
        RETURN
      ENDIF
C
C---- Newton iteration to get exact SLE value
      DO 20 ITER=1, 50
        XLE  = SEVAL(SLE,X,XP,S,N)
        YLE  = SEVAL(SLE,Y,YP,S,N)
        DXDS = DEVAL(SLE,X,XP,S,N)
        DYDS = DEVAL(SLE,Y,YP,S,N)
        DXDD = D2VAL(SLE,X,XP,S,N)
        DYDD = D2VAL(SLE,Y,YP,S,N)
C
        XCHORD = XLE - XTE
        YCHORD = YLE - YTE
C
C------ drive dot product between chord line and LE tangent to zero
        RES  = XCHORD*DXDS + YCHORD*DYDS
        RESS = DXDS  *DXDS + DYDS  *DYDS
     &       + XCHORD*DXDD + YCHORD*DYDD
C
C------ Newton delta for SLE 
        DSLE = -RES/RESS
C
        DSLE = MAX( DSLE , -0.02*ABS(XCHORD+YCHORD) )
        DSLE = MIN( DSLE ,  0.02*ABS(XCHORD+YCHORD) )
        SLE = SLE + DSLE
        IF(ABS(DSLE) .LT. DSEPS) RETURN
   20 CONTINUE
C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE)
     &  WRITE(*,*) 'LEFIND:  LE point not found.  Continuing...'
      SLE = S(I)
      RETURN
      END

C===================================================================70
C
C     Calculates curvature of splined 2-D curve |
C     at S = SS                                 |
C                                               |
C     S        arc length array of curve        |
C     X, Y     coordinate arrays of curve       |
C     XS,YS    derivative arrays                |
C              (calculated earlier by SPLINE)   |
C
C===================================================================70
      FUNCTION CURV(SS,X,XS,Y,YS,S,N)
      DIMENSION X(N), XS(N), Y(N), YS(N), S(N)
C     
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
C
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      XD = X(I) - X(I-1) + (1.0-4.0*T+3.0*T*T)*CX1 + T*(3.0*T-2.0)*CX2
      XDD = (6.0*T-4.0)*CX1 + (6.0*T-2.0)*CX2
C
      CY1 = DS*YS(I-1) - Y(I) + Y(I-1)
      CY2 = DS*YS(I)   - Y(I) + Y(I-1)
      YD = Y(I) - Y(I-1) + (1.0-4.0*T+3.0*T*T)*CY1 + T*(3.0*T-2.0)*CY2
      YDD = (6.0*T-4.0)*CY1 + (6.0*T-2.0)*CY2
C 
      SD = SQRT(XD*XD + YD*YD)
      SD = MAX(SD,0.001*DS)
C
      CURV = (XD*YDD - YD*XDD) / SD**3
C
      RETURN
      END ! CURV

C===================================================================70
C
C     Calculates X(SS)                             |
C     XS array must have been calculated by SPLINE |
C
C===================================================================70
      FUNCTION SEVAL(SS,X,XS,S,N)
      DIMENSION X(N), XS(N), S(N)
      ILOW = 1
      I = N
C
   10 IF(I-ILOW .LE. 1) GO TO 11
C
      IMID = (I+ILOW)/2
      IF(SS .LT. S(IMID)) THEN
       I = IMID
      ELSE
       ILOW = IMID
      ENDIF
      GO TO 10
C
   11 DS = S(I) - S(I-1)
      T = (SS - S(I-1)) / DS
      CX1 = DS*XS(I-1) - X(I) + X(I-1)
      CX2 = DS*XS(I)   - X(I) + X(I-1)
      SEVAL = T*X(I) + (1.0-T)*X(I-1) + (T-T*T)*((1.0-T)*CX1 - T*CX2)
      RETURN
      END ! SEVAL

C===================================================================70
C
C     Calculates normal unit vector
C     components at airfoil panel nodes
C
C===================================================================70
      SUBROUTINE NCALC(X,Y,S,N,XN,YN)
      DIMENSION X(N), Y(N), S(N), XN(N), YN(N)
C
      IF(N.LE.1) RETURN
C
      CALL SEGSPL(X,XN,S,N)
      CALL SEGSPL(Y,YN,S,N)
      DO 10 I=1, N
        SX =  YN(I)
        SY = -XN(I)
        SMOD = SQRT(SX*SX + SY*SY)
        XN(I) = SX/SMOD
        YN(I) = SY/SMOD
   10 CONTINUE
C
C---- average normal vectors at corner points
      DO 20 I=1, N-1
        IF(S(I) .EQ. S(I+1)) THEN
          SX = 0.5*(XN(I) + XN(I+1))
          SY = 0.5*(YN(I) + YN(I+1))
          SMOD = SQRT(SX*SX + SY*SY)
          XN(I)   = SX/SMOD
          YN(I)   = SY/SMOD
          XN(I+1) = SX/SMOD
          YN(I+1) = SY/SMOD
        ENDIF
 20   CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     set angles of airfoil panels
C
C===================================================================70
      SUBROUTINE APCALC

      use xfoil_inc
C
      DO 10 I=1, N-1
        SX = X(I+1) - X(I)
        SY = Y(I+1) - Y(I)
        IF(SX.EQ.0.0 .AND. SY.EQ.0.0) THEN
          APANEL(I) = ATAN2( -NY(I) , -NX(I) )
        ELSE
          APANEL(I) = ATAN2( SX , -SY )
        ENDIF
   10 CONTINUE
C
C---- TE panel
      I = N
      IP = 1
      IF(SHARP) THEN
       APANEL(I) = PI
      ELSE
       SX = X(IP) - X(I)
       SY = Y(IP) - Y(I)
       APANEL(I) = ATAN2( -SX , SY ) + PI
      ENDIF
C
      RETURN
      END

C===================================================================70
C
C     Calculates arc length SOPP of point 
C     which is opposite of point SI, on the 
C     other side of the airfoil baseline
C
C===================================================================70
      SUBROUTINE SOPPS(SOPP, SI, X,XP,Y,YP,S,N, SLE, SILENT_MODE)

      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
      LOGICAL :: SILENT_MODE
C
C---- reference length for testing convergence
      SLEN = S(N) - S(1)
C
C---- set chordline vector
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHORD = SQRT((XTE-XLE)**2 + (YTE-YLE)**2)
      DXC = (XTE-XLE) / CHORD
      DYC = (YTE-YLE) / CHORD
C
      IF(SI.LT.SLE) THEN
       IN = 1
       INOPP = N
      ELSE
       IN = N
       INOPP = 1
      ENDIF
      SFRAC = (SI-SLE)/(S(IN)-SLE)
      SOPP = SLE + SFRAC*(S(INOPP)-SLE)
C     
      IF(ABS(SFRAC) .LE. 1.0E-5) THEN
       SOPP = SLE
       RETURN
      ENDIF
C
C---- XBAR = x coordinate in chord-line axes
      XI  = SEVAL(SI , X,XP,S,N)
      YI  = SEVAL(SI , Y,YP,S,N)
      XLE = SEVAL(SLE, X,XP,S,N)
      YLE = SEVAL(SLE, Y,YP,S,N)
      XBAR = (XI-XLE)*DXC + (YI-YLE)*DYC
C
C---- converge on exact opposite point with same XBAR value
      DO 300 ITER=1, 12
        XOPP  = SEVAL(SOPP,X,XP,S,N)
        YOPP  = SEVAL(SOPP,Y,YP,S,N)
        XOPPD = DEVAL(SOPP,X,XP,S,N)
        YOPPD = DEVAL(SOPP,Y,YP,S,N)
C
        RES  = (XOPP -XLE)*DXC + (YOPP -YLE)*DYC - XBAR
        RESD =  XOPPD     *DXC +  YOPPD     *DYC
C
        IF(ABS(RES)/SLEN .LT. 1.0E-5) GO TO 305
        IF(RESD .EQ. 0.0) GO TO 303
C
        DSOPP = -RES/RESD
        SOPP = SOPP + DSOPP
C
        IF(ABS(DSOPP)/SLEN .LT. 1.0E-5) GO TO 305
 300  CONTINUE
C     DP mod: added SILENT_MODE option
 303  IF (.NOT. SILENT_MODE) WRITE(*,*)
     &      'SOPPS: Opposite-point location failed. Continuing...'
      SOPP = SLE + SFRAC*(S(INOPP)-SLE)
C
 305  CONTINUE
      RETURN
      END ! SOPPS

C===================================================================70
C
C     Calculates max thickness and camber at airfoil points
C     DP mod: calculates min thickness aft of x = 0.5
C
C     Note: this routine does not find the maximum camber or 
C           thickness exactly as it only looks at discrete points
C
C     Input:
C       N      number of points
C       X(.)   shape coordinate point arrays
C       Y(.)
C
C     Output:
C       THICK  max thickness
C       CAMBR  max camber
C       THICKM min thickness aft of x = 0.5
C
C===================================================================70
      SUBROUTINE TCCALC(X,XP,Y,YP,S,N,SILENT_MODE,
     &                  THICK,XTHICK,THICKM,XTHICKM,CAMBR,XCAMBR)

      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
      LOGICAL :: SILENT_MODE
      CALL LEFIND(SLE,X,XP,Y,YP,S,N,SILENT_MODE)
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHORD = SQRT((XTE-XLE)**2 + (YTE-YLE)**2)
C
C---- set unit chord-line vector
      DXC = (XTE-XLE) / CHORD
      DYC = (YTE-YLE) / CHORD
C
      THICK = 0.
      XTHICK = 0.
      THICKM = 1000.
      XTHICKM = 0.
      CAMBR = 0.
      XCAMBR = 0.
C
C---- go over each point, finding the y-thickness and camber
      DO 30 I=1, N
        XBAR = (X(I)-XLE)*DXC + (Y(I)-YLE)*DYC
        YBAR = (Y(I)-YLE)*DXC - (X(I)-XLE)*DYC
C
C------ set point on the opposite side with the same chord x value
        CALL SOPPS(SOPP, S(I), X,XP,Y,YP,S,N, SLE, SILENT_MODE)
        XOPP = SEVAL(SOPP,X,XP,S,N)
        YOPP = SEVAL(SOPP,Y,YP,S,N)
C
        YBAROP = (YOPP-YLE)*DXC - (XOPP-XLE)*DYC
C
        YC = 0.5*(YBAR+YBAROP)
        YT =  ABS(YBAR-YBAROP)
C
        IF(ABS(YC) .GT. ABS(CAMBR)) THEN
         CAMBR = YC
         XCAMBR = XOPP
        ENDIF
        IF(ABS(YT) .GT. ABS(THICK)) THEN
         THICK = YT
         XTHICK = XOPP
        ENDIF
C       DP mod: calculates min thickness aft of x = 0.5
        IF((ABS(YT) .LT. ABS(THICKM)) .AND. (XOPP > 0.5)) THEN
          THICKM = YT
          XTHICKM = XOPP
        ENDIF
   30 CONTINUE
C
      RETURN
      END ! TCCALC

C===================================================================70
C
C     IPRINT=2:   Displays all panel node corner angles
C     IPRINT=1:   Displays max panel node corner angle
C     IPRINT=0:   No display... just returns values
C
C===================================================================70
      SUBROUTINE CANG(X,Y,N,IPRINT, IMAX,AMAX)
      DIMENSION X(*), Y(*)
C
      AMAX = 0.0
      IMAX = 1
C
C---- go over each point, calculating corner angle
      IF(IPRINT.EQ.2) WRITE(*,1050)
      DO 30 I=2, N-1
        DX1 = X(I) - X(I-1)
        DY1 = Y(I) - Y(I-1)
        DX2 = X(I) - X(I+1)
        DY2 = Y(I) - Y(I+1)
C
C------ allow for doubled points
        IF(DX1.EQ.0.0 .AND. DY1.EQ.0.0) THEN
         DX1 = X(I) - X(I-2)
         DY1 = Y(I) - Y(I-2)
        ENDIF
        IF(DX2.EQ.0.0 .AND. DY2.EQ.0.0) THEN
         DX2 = X(I) - X(I+2)
         DY2 = Y(I) - Y(I+2)
        ENDIF
C
        CROSSP = (DX2*DY1 - DY2*DX1)
     &         / SQRT((DX1**2 + DY1**2) * (DX2**2 + DY2**2))
        ANGL = ASIN(CROSSP)*(180.0/3.1415926)
        IF(IPRINT.EQ.2) WRITE(*,1100) I, X(I), Y(I), ANGL
        IF(ABS(ANGL) .GT. ABS(AMAX)) THEN
         AMAX = ANGL
         IMAX = I
        ENDIF
   30 CONTINUE
C
      IF(IPRINT.GE.1) WRITE(*,1200) AMAX, IMAX, X(IMAX), Y(IMAX)
C
      RETURN
C
 1050 FORMAT(/'  i       x        y      angle')
CCC             120   0.2134  -0.0234   25.322
 1100 FORMAT(1X,I3, 2F9.4, F9.3)
 1200 FORMAT(/' Maximum panel corner angle =', F7.3,
     &        '   at  i,x,y  = ', I3, 2F9.4 )
      END ! CANG

C===================================================================70
C
C     Set paneling distribution from buffer airfoil
C     geometry, thus creating current airfoil.
C 
C     If REFINE=True, bunch points at x=XSREF on
C     top side and at x=XPREF on bottom side
C     by setting a fictitious local curvature of
C     CTRRAT*(LE curvature) there.
C
C===================================================================70
      SUBROUTINE PANGEN(SHOPAR)

      use xfoil_inc

      LOGICAL SHOPAR
C
C     DP mod: added SILENT_MODE option
      IF(NB.LT.2) THEN
       IF (.NOT. SILENT_MODE)
     &   WRITE(*,*) 'PANGEN: Buffer airfoil not available.'
       N = 0
       RETURN
      ENDIF
C
C---- Number of temporary nodes for panel distribution calculation
C       exceeds the specified panel number by factor of IPFAC.
      IPFAC = 3
      IPFAC = 5
C
C---- number of airfoil panel points
      N = NPAN
C
cC---- number of wake points
c      NW = NPAN/8 + 2
c      IF(NW.GT.IWX) THEN
c       WRITE(*,*)
c     &  'Array size (IWX) too small.  Last wake point index reduced.'
c       NW = IWX
c      ENDIF
C
C---- set arc length spline parameter
      CALL SCALC(XB,YB,SB,NB)
C
C---- spline raw airfoil coordinates
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
C
C---- normalizing length (~ chord)
      SBREF = 0.5*(SB(NB)-SB(1))
C
C---- set up curvature array
      DO I = 1, NB
        W5(I) = ABS( CURV(SB(I),XB,XBP,YB,YBP,SB,NB) ) * SBREF
      ENDDO
C
C---- locate LE point arc length value and the normalized curvature there
      CALL LEFIND(SBLE,XB,XBP,YB,YBP,SB,NB,SILENT_MODE)
      CVLE = ABS( CURV(SBLE,XB,XBP,YB,YBP,SB,NB) ) * SBREF
C
C---- check for doubled point (sharp corner) at LE
      IBLE = 0
      DO I = 1, NB-1
        IF(SBLE.EQ.SB(I) .AND. SBLE.EQ.SB(I+1)) THEN
         IBLE = I
C        DP mod: added SILENT_MODE option
         IF (.NOT. SILENT_MODE) THEN
           WRITE(*,*)
           WRITE(*,*) 'Sharp leading edge'
         ENDIF
         GO TO 21
        ENDIF
      ENDDO
 21   CONTINUE
C
C---- set LE, TE points
      XBLE = SEVAL(SBLE,XB,XBP,SB,NB)
      YBLE = SEVAL(SBLE,YB,YBP,SB,NB)
      XBTE = 0.5*(XB(1)+XB(NB))
      YBTE = 0.5*(YB(1)+YB(NB))
      CHBSQ = (XBTE-XBLE)**2 + (YBTE-YBLE)**2
C
C---- set average curvature over 2*NK+1 points within Rcurv of LE point
      NK = 3
      CVSUM = 0.
      DO K = -NK, NK
        FRAC = FLOAT(K)/FLOAT(NK)
        SBK = SBLE + FRAC*SBREF/MAX(CVLE,20.0)
        CVK = ABS( CURV(SBK,XB,XBP,YB,YBP,SB,NB) ) * SBREF
        CVSUM = CVSUM + CVK
      ENDDO
      CVAVG = CVSUM/FLOAT(2*NK+1)
C
C---- dummy curvature for sharp LE
      IF(IBLE.NE.0) CVAVG = 10.0
C
C---- set curvature attraction coefficient actually used
      CC = 6.0 * CVPAR
C
C---- set artificial curvature at TE to bunch panels there
      CVTE = CVAVG * CTERAT
      W5(1)  = CVTE
      W5(NB) = CVTE
C
C
C**** smooth curvature array for smoother panel size distribution  ****
C
CCC      CALL ASKR('Enter curvature smoothing length/c^',SMOOL)
CCC      SMOOL = 0.010
C
C---- set smoothing length = 1 / averaged LE curvature, but 
C-    no more than 5% of chord and no less than 1/4 average panel spacing
      SMOOL = MAX( 1.0/MAX(CVAVG,20.0) , 0.25 /FLOAT(NPAN/2) )
C
      SMOOSQ = (SMOOL*SBREF) ** 2
C
C---- set up tri-diagonal system for smoothed curvatures
      W2(1) = 1.0
      W3(1) = 0.0
      DO I=2, NB-1
        DSM = SB(I) - SB(I-1)
        DSP = SB(I+1) - SB(I)
        DSO = 0.5*(SB(I+1) - SB(I-1))
C
        IF(DSM.EQ.0.0 .OR. DSP.EQ.0.0) THEN
C------- leave curvature at corner point unchanged
         W1(I) = 0.0
         W2(I) = 1.0
         W3(I) = 0.0
        ELSE
         W1(I) =  SMOOSQ * (         - 1.0/DSM) / DSO
         W2(I) =  SMOOSQ * ( 1.0/DSP + 1.0/DSM) / DSO  +  1.0
         W3(I) =  SMOOSQ * (-1.0/DSP          ) / DSO
        ENDIF
      ENDDO
C
      W1(NB) = 0.0
      W2(NB) = 1.0
C
C---- fix curvature at LE point by modifying equations adjacent to LE
      DO I=2, NB-1
        IF(SB(I).EQ.SBLE .OR. I.EQ.IBLE .OR. I.EQ.IBLE+1) THEN
C------- if node falls right on LE point, fix curvature there
         W1(I) = 0.
         W2(I) = 1.0
         W3(I) = 0.
         W5(I) = CVLE
        ELSE IF(SB(I-1).LT.SBLE .AND. SB(I).GT.SBLE) THEN
C------- modify equation at node just before LE point
         DSM = SB(I-1) - SB(I-2)
         DSP = SBLE    - SB(I-1)
         DSO = 0.5*(SBLE - SB(I-2))
C
         W1(I-1) =  SMOOSQ * (         - 1.0/DSM) / DSO
         W2(I-1) =  SMOOSQ * ( 1.0/DSP + 1.0/DSM) / DSO  +  1.0
         W3(I-1) =  0.
         W5(I-1) = W5(I-1) + SMOOSQ*CVLE/(DSP*DSO)
C
C------- modify equation at node just after LE point
         DSM = SB(I) - SBLE
         DSP = SB(I+1) - SB(I)
         DSO = 0.5*(SB(I+1) - SBLE)
         W1(I) =  0.
         W2(I) =  SMOOSQ * ( 1.0/DSP + 1.0/DSM) / DSO  +  1.0
         W3(I) =  SMOOSQ * (-1.0/DSP          ) / DSO
         W5(I) = W5(I) + SMOOSQ*CVLE/(DSM*DSO)
C
         GO TO 51
        ENDIF
      ENDDO
   51 CONTINUE
C
C---- set artificial curvature at bunching points and fix it there
      DO I=2, NB-1
C------ chord-based x/c coordinate
        XOC = (  (XB(I)-XBLE)*(XBTE-XBLE)
     &         + (YB(I)-YBLE)*(YBTE-YBLE) ) / CHBSQ
C
        IF(SB(I).LT.SBLE) THEN
C------- check if top side point is in refinement area
         IF(XOC.GT.XSREF1 .AND. XOC.LT.XSREF2) THEN
          W1(I) = 0.
          W2(I) = 1.0
          W3(I) = 0.
          W5(I) = CVLE*CTRRAT
         ENDIF
        ELSE
C------- check if bottom side point is in refinement area
         IF(XOC.GT.XPREF1 .AND. XOC.LT.XPREF2) THEN
          W1(I) = 0.
          W2(I) = 1.0
          W3(I) = 0.
          W5(I) = CVLE*CTRRAT
         ENDIF
        ENDIF
      ENDDO
C
C---- solve for smoothed curvature array W5
      IF(IBLE.EQ.0) THEN
       CALL TRISOL(W2,W1,W3,W5,NB)
      ELSE
       I = 1
       CALL TRISOL(W2(I),W1(I),W3(I),W5(I),IBLE)
       I = IBLE+1
       CALL TRISOL(W2(I),W1(I),W3(I),W5(I),NB-IBLE)
      ENDIF
C
C---- find max curvature
      CVMAX = 0.
      DO I=1, NB
        CVMAX = MAX( CVMAX , ABS(W5(I)) )
      ENDDO
C
C---- normalize curvature array
      DO I=1, NB
        W5(I) = W5(I) / CVMAX
      ENDDO
C
C---- spline curvature array
      CALL SEGSPL(W5,W6,SB,NB)
C
C---- Set initial guess for node positions uniform in s.
C     More nodes than specified (by factor of IPFAC) are 
C     temporarily used  for more reliable convergence.
      NN = IPFAC*(N-1)+1
C
C---- ratio of lengths of panel at TE to one away from the TE
      RDSTE = 0.667
      RTF = (RDSTE-1.0)*2.0 + 1.0
C
      IF(IBLE.EQ.0) THEN
C
       DSAVG = (SB(NB)-SB(1))/(FLOAT(NN-3) + 2.0*RTF)
       SNEW(1) = SB(1)
       DO I=2, NN-1
         SNEW(I) = SB(1) + DSAVG * (FLOAT(I-2) + RTF)
       ENDDO
       SNEW(NN) = SB(NB)
C
      ELSE
C
       NFRAC1 = (N * IBLE) / NB
C
       NN1 = IPFAC*(NFRAC1-1)+1
       DSAVG1 = (SBLE-SB(1))/(FLOAT(NN1-2) + RTF)
       SNEW(1) = SB(1)
       DO I=2, NN1
         SNEW(I) = SB(1) + DSAVG1 * (FLOAT(I-2) + RTF)
       ENDDO
C
       NN2 = NN - NN1 + 1
       DSAVG2 = (SB(NB)-SBLE)/(FLOAT(NN2-2) + RTF)
       DO I=2, NN2-1
         SNEW(I-1+NN1) = SBLE + DSAVG2 * (FLOAT(I-2) + RTF)
       ENDDO
       SNEW(NN) = SB(NB)
C
      ENDIF
C
C---- Newton iteration loop for new node positions
      DO 10 ITER=1, 20
C
C------ set up tri-diagonal system for node position deltas
        CV1  = SEVAL(SNEW(1),W5,W6,SB,NB)
        CV2  = SEVAL(SNEW(2),W5,W6,SB,NB)
        CVS1 = DEVAL(SNEW(1),W5,W6,SB,NB)
        CVS2 = DEVAL(SNEW(2),W5,W6,SB,NB)
C
        CAVM = SQRT(CV1**2 + CV2**2)
        IF(CAVM .EQ. 0.0) THEN
          CAVM_S1 = 0.
          CAVM_S2 = 0.
        ELSE
          CAVM_S1 = CVS1 * CV1/CAVM
          CAVM_S2 = CVS2 * CV2/CAVM
        ENDIF
C
        DO 110 I=2, NN-1
          DSM = SNEW(I) - SNEW(I-1)
          DSP = SNEW(I) - SNEW(I+1)
          CV3  = SEVAL(SNEW(I+1),W5,W6,SB,NB)
          CVS3 = DEVAL(SNEW(I+1),W5,W6,SB,NB)
C
          CAVP = SQRT(CV3**2 + CV2**2)
          IF(CAVP .EQ. 0.0) THEN
            CAVP_S2 = 0.
            CAVP_S3 = 0.
          ELSE
            CAVP_S2 = CVS2 * CV2/CAVP
            CAVP_S3 = CVS3 * CV3/CAVP
          ENDIF
C
          FM = CC*CAVM + 1.0
          FP = CC*CAVP + 1.0
C
          REZ = DSP*FP + DSM*FM
C
C-------- lower, main, and upper diagonals
          W1(I) =      -FM  +  CC*               DSM*CAVM_S1
          W2(I) =  FP + FM  +  CC*(DSP*CAVP_S2 + DSM*CAVM_S2)
          W3(I) = -FP       +  CC* DSP*CAVP_S3
C
C-------- residual, requiring that
C         (1 + C*curv)*deltaS is equal on both sides of node i
          W4(I) = -REZ
C
          CV1 = CV2
          CV2 = CV3
          CVS1 = CVS2
          CVS2 = CVS3
          CAVM    = CAVP
          CAVM_S1 = CAVP_S2
          CAVM_S2 = CAVP_S3
  110   CONTINUE
C
C------ fix endpoints (at TE)
        W2(1) = 1.0
        W3(1) = 0.0
        W4(1) = 0.0
        W1(NN) = 0.0
        W2(NN) = 1.0
        W4(NN) = 0.0
C
        IF(RTF .NE. 1.0) THEN
C------- fudge equations adjacent to TE to get TE panel length ratio RTF
C
         I = 2
         W4(I) = -((SNEW(I) - SNEW(I-1)) + RTF*(SNEW(I) - SNEW(I+1)))
         W1(I) = -1.0
         W2(I) =  1.0 + RTF
         W3(I) =      - RTF
C
         I = NN-1
         W4(I) = -((SNEW(I) - SNEW(I+1)) + RTF*(SNEW(I) - SNEW(I-1)))
         W3(I) = -1.0
         W2(I) =  1.0 + RTF
         W1(I) =      - RTF
        ENDIF
C
C
C------ fix sharp LE point
        IF(IBLE.NE.0) THEN
         I = NN1
         W1(I) = 0.0
         W2(I) = 1.0
         W3(I) = 0.0
         W4(I) = SBLE - SNEW(I)
        ENDIF
C
C------ solve for changes W4 in node position arc length values
        CALL TRISOL(W2,W1,W3,W4,NN)
C
C------ find under-relaxation factor to keep nodes from changing order
        RLX = 1.0
        DMAX = 0.0
        DO I=1, NN-1
          DS  = SNEW(I+1) - SNEW(I)
          DDS = W4(I+1) - W4(I)
          DSRAT = 1.0 + RLX*DDS/DS
          IF(DSRAT.GT.4.0) RLX = (4.0-1.0)*DS/DDS
          IF(DSRAT.LT.0.2) RLX = (0.2-1.0)*DS/DDS
          DMAX = MAX(ABS(W4(I)),DMAX)
        ENDDO
C
C------ update node position
        DO I=2, NN-1
          SNEW(I) = SNEW(I) + RLX*W4(I)
        ENDDO
C
CCC        IF(RLX.EQ.1.0) WRITE(*,*) DMAX
CCC        IF(RLX.NE.1.0) WRITE(*,*) DMAX,'    RLX =',RLX
        IF(ABS(DMAX).LT.1.E-3) GO TO 11
   10 CONTINUE
C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE) 
     &  WRITE(*,*) 'Paneling convergence failed.  Continuing anyway...'
C
   11 CONTINUE
C
C---- set new panel node coordinates
      DO I=1, N
        IND = IPFAC*(I-1) + 1
        S(I) = SNEW(IND)
        X(I) = SEVAL(SNEW(IND),XB,XBP,SB,NB)
        Y(I) = SEVAL(SNEW(IND),YB,YBP,SB,NB)
      ENDDO
C
C
C---- go over buffer airfoil again, checking for corners (double points)
      NCORN = 0
      DO 25 IB=1, NB-1
        IF(SB(IB) .EQ. SB(IB+1)) THEN
C------- found one !
C
         NCORN = NCORN+1
         XBCORN = XB(IB)
         YBCORN = YB(IB)
         SBCORN = SB(IB)
C
C------- find current-airfoil panel which contains corner
         DO 252 I=1, N
C
C--------- keep stepping until first node past corner
           IF(S(I) .LE. SBCORN) GO TO 252
C
C---------- move remainder of panel nodes to make room for additional node
            DO 2522 J=N, I, -1
              X(J+1) = X(J)
              Y(J+1) = Y(J)
              S(J+1) = S(J)
 2522       CONTINUE
            N = N+1
C
            IF(N .GT. IQX-1)
     &       STOP 'PANEL: Too many panels. Increase IQX in xfoil_inc'
C
            X(I) = XBCORN
            Y(I) = YBCORN
            S(I) = SBCORN
C
C---------- shift nodes adjacent to corner to keep panel sizes comparable
            IF(I-2 .GE. 1) THEN
             S(I-1) = 0.5*(S(I) + S(I-2))
             X(I-1) = SEVAL(S(I-1),XB,XBP,SB,NB)
             Y(I-1) = SEVAL(S(I-1),YB,YBP,SB,NB)
            ENDIF
C
            IF(I+2 .LE. N) THEN
             S(I+1) = 0.5*(S(I) + S(I+2))
             X(I+1) = SEVAL(S(I+1),XB,XBP,SB,NB)
             Y(I+1) = SEVAL(S(I+1),YB,YBP,SB,NB)
            ENDIF
C
C---------- go on to next input geometry point to check for corner
            GO TO 25
C
  252    CONTINUE
        ENDIF
   25 CONTINUE
C
      CALL SCALC(X,Y,S,N)
      CALL SEGSPL(X,XP,S,N)
      CALL SEGSPL(Y,YP,S,N)
      CALL LEFIND(SLE,X,XP,Y,YP,S,N,SILENT_MODE)
C
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHORD  = SQRT( (XTE-XLE)**2 + (YTE-YLE)**2 )
C
C---- calculate panel size ratios (user info)
      DSMIN =  1000.0
      DSMAX = -1000.0
      DO 40 I=1, N-1
        DS = S(I+1)-S(I)
        IF(DS .EQ. 0.0) GO TO 40
          DSMIN = MIN(DSMIN,DS)
          DSMAX = MAX(DSMAX,DS)
   40 CONTINUE
C
      DSMIN = DSMIN*FLOAT(N-1)/S(N)
      DSMAX = DSMAX*FLOAT(N-1)/S(N)
ccc      WRITE(*,*) 'DSmin/DSavg = ',DSMIN,'     DSmax/DSavg = ',DSMAX
C
C---- set various flags for new airfoil
      LGAMU = .FALSE.
      LWAKE = .FALSE.
      LQAIJ = .FALSE.
      LADIJ = .FALSE.
      LWDIJ = .FALSE.
      LIPAN = .FALSE.
      LBLINI = .FALSE.
      LVCONV = .FALSE.
C
C     DP note: commented these lines (not using flap)
C      IF(LBFLAP) THEN
C       XOF = XBF
C       YOF = YBF
C       LFLAP = .TRUE.
C      ENDIF
C
C---- determine if TE is blunt or sharp, calculate TE geometry parameters
      CALL TECALC
C
C---- calculate normal vectors
      CALL NCALC(X,Y,S,N,NX,NY)
C
C---- calculate panel angles for panel routines
      CALL APCALC
C
C     DP mod: added SILENT_MODE option
      IF(SHARP) THEN
       IF (.NOT. SILENT_MODE) WRITE(*,1090) 'Sharp trailing edge'
      ELSE
       GAP = SQRT((X(1)-X(N))**2 + (Y(1)-Y(N))**2)
       IF (.NOT. SILENT_MODE) 
     &   WRITE(*,1090) 'Blunt trailing edge.  Gap =', GAP
      ENDIF
 1090 FORMAT(/1X,A,F9.5)
C
C     DP mod: added SILENT_MODE option
      IF(SHOPAR .AND. .NOT. SILENT_MODE) 
     &           WRITE(*,1100) NPAN, CVPAR, CTERAT, CTRRAT,
     &                         XSREF1, XSREF2, XPREF1, XPREF2
 1100 FORMAT(/' Paneling parameters used...'
     &       /'   Number of panel nodes      ' , I4
     &       /'   Panel bunching parameter   ' , F6.3
     &       /'   TE/LE panel density ratio  ' , F6.3
     &       /'   Refined-area/LE panel density ratio   ' , F6.3
     &       /'   Top    side refined area x/c limits ' , 2F6.3
     &       /'   Bottom side refined area x/c limits ' , 2F6.3)

C     DP mod: added thickness and camber calculations here
      CALL TCCALC(X,XP,Y,YP,S,N,SILENT_MODE,
     &            THICKB,XTHICKB,THICKM,XTHICKM,CAMBR,XCAMBR)
C
C     DP mod: added panel corner angle calculations here
      CALL CANG(X,Y,N,0,IMAX,AMAX)
C
      RETURN
      END ! PANGEN

C===================================================================70
C
C     Sets geometrically stretched array S:
C
C       S(i+1) - S(i)  =  r * [S(i) - S(i-1)]
C
C       S     (output)  array to be set  
C       DS1   (input)   first S increment:  S(2) - S(1)
C       SMAX  (input)   final S value:      S(NN)
C       NN    (input)   number of points
C
C===================================================================70
      SUBROUTINE SETEXP(S,DS1,SMAX,NN,SILENT_MODE)

      REAL*8 S(NN)
      LOGICAL SILENT_MODE
C
      SIGMA = SMAX/DS1
      NEX = NN-1
      RNEX = FLOAT(NEX)
      RNI = 1.0/RNEX
C
C---- solve quadratic for initial geometric ratio guess
      AAA = RNEX*(RNEX-1.0)*(RNEX-2.0) / 6.0
      BBB = RNEX*(RNEX-1.0) / 2.0
      CCC = RNEX - SIGMA
C
      DISC = BBB**2 - 4.0*AAA*CCC
      DISC = MAX( 0.0 , DISC )
C
      IF(NEX.LE.1) THEN
       STOP 'SETEXP: Cannot fill array.  N too small.'
      ELSE IF(NEX.EQ.2) THEN
       RATIO = -CCC/BBB  +  1.0
      ELSE
       RATIO = (-BBB + SQRT(DISC))/(2.0*AAA)  +  1.0
      ENDIF
C
      IF(RATIO.EQ.1.0) GO TO 11
C
C---- Newton iteration for actual geometric ratio
      DO 1 ITER=1, 100
        SIGMAN = (RATIO**NEX - 1.0) / (RATIO - 1.0)
        RES = SIGMAN**RNI - SIGMA**RNI
        DRESDR = RNI*SIGMAN**RNI
     &         * (RNEX*RATIO**(NEX-1) - SIGMAN) / (RATIO**NEX - 1.0)
C
        DRATIO = -RES/DRESDR
        RATIO = RATIO + DRATIO
C
        IF(ABS(DRATIO) .LT. 1.0E-5) GO TO 11
C
    1 CONTINUE
C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE)
     &  WRITE(*,*) 'SETEXP: Convergence failed.  Continuing anyway ...'
C
C---- set up stretched array using converged geometric ratio
   11 S(1) = 0.0
      DS = DS1
      DO 2 N=2, NN
        S(N) = S(N-1) + DS
        DS = DS*RATIO
    2 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Sets wake coordinate array for current surface 
C     vorticity and/or mass source distributions.
C
C===================================================================70
      SUBROUTINE XYWAKE

      use xfoil_inc
C
C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE)
     &  WRITE(*,*) 'Calculating wake trajectory ...'
C
C---- number of wake points
      NW = N/8 + 2
      IF(NW.GT.IWX) THEN
C      DP mod: added SILENT_MODE option
       IF (.NOT. SILENT_MODE) WRITE(*,*)
     &  'Array size (IWX) too small.  Last wake point index reduced.'
       NW = IWX
      ENDIF
C
      DS1 = 0.5*(S(2) - S(1) + S(N) - S(N-1))
      CALL SETEXP(SNEW(N+1),DS1,WAKLEN*CHORD,NW,SILENT_MODE)
C
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
C
C---- set first wake point a tiny distance behind TE
      I = N+1
      SX = 0.5*(YP(N) - YP(1))
      SY = 0.5*(XP(1) - XP(N))
      SMOD = SQRT(SX**2 + SY**2)
      NX(I) = SX / SMOD
      NY(I) = SY / SMOD
      X(I) = XTE - 0.0001*NY(I)
      Y(I) = YTE + 0.0001*NX(I)
      S(I) = S(N)
C
C---- calculate streamfunction gradient components at first point
      CALL PSILIN(I,X(I),Y(I),1.0,0.0,PSI,PSI_X,.FALSE.,.FALSE.)
      CALL PSILIN(I,X(I),Y(I),0.0,1.0,PSI,PSI_Y,.FALSE.,.FALSE.)
C
C---- set unit vector normal to wake at first point
      NX(I+1) = -PSI_X / SQRT(PSI_X**2 + PSI_Y**2)
      NY(I+1) = -PSI_Y / SQRT(PSI_X**2 + PSI_Y**2)
C
C---- set angle of wake panel normal
      APANEL(I) = ATAN2( PSI_Y , PSI_X )
C
C---- set rest of wake points
      DO 10 I=N+2, N+NW
        DS = SNEW(I) - SNEW(I-1)
C
C------ set new point DS downstream of last point
        X(I) = X(I-1) - DS*NY(I)
        Y(I) = Y(I-1) + DS*NX(I)
        S(I) = S(I-1) + DS
C
        IF(I.EQ.N+NW) GO TO 10
C
C------- calculate normal vector for next point
         CALL PSILIN(I,X(I),Y(I),1.0,0.0,PSI,PSI_X,.FALSE.,.FALSE.)
         CALL PSILIN(I,X(I),Y(I),0.0,1.0,PSI,PSI_Y,.FALSE.,.FALSE.)
C
         NX(I+1) = -PSI_X / SQRT(PSI_X**2 + PSI_Y**2)
         NY(I+1) = -PSI_Y / SQRT(PSI_X**2 + PSI_Y**2)
C
C------- set angle of wake panel normal
         APANEL(I) = ATAN2( PSI_Y , PSI_X )
C
   10 CONTINUE
C
C---- set wake presence flag and corresponding alpha
      LWAKE = .TRUE.
      AWAKE =  ALFA
C
C---- old source influence matrix is invalid for the new wake geometry
      LWDIJ = .FALSE.
C
      RETURN
      END

C===================================================================70
C
C     Calculates current streamfunction Psi and tangential velocity
C     Qtan at panel node or wake node I due to freestream and wake
C     sources Sig.  Also calculates sensitivity vectors dPsi/dSig
C     (DZDM) and dQtan/dSig (DQDM).
C
C          Airfoil:  1   < I < N
C          Wake:     N+1 < I < N+NW
C
C===================================================================70
      SUBROUTINE PSWLIN(I,XI,YI,NXI,NYI,PSI,PSI_NI)

      use xfoil_inc

      REAL*8 :: NXI, NYI
C
      IO = I
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
C
      DO 4 JO=N+1, N+NW
        DZDM(JO) = 0.0
        DQDM(JO) = 0.0
    4 CONTINUE
C
      PSI    = 0.
      PSI_NI = 0.
C
      DO 20 JO=N+1, N+NW-1
C
        JP = JO+1
C
        JM = JO-1
        JQ = JP+1
        IF(JO.EQ.N+1) THEN
         JM = JO
        ELSE IF(JO.EQ.N+NW-1) THEN
         JQ = JP
        ENDIF
C
        DSO = SQRT((X(JO)-X(JP))**2 + (Y(JO)-Y(JP))**2)
        DSIO = 1.0 / DSO
C
        APAN = APANEL(JO)
C
        RX1 = XI - X(JO)
        RY1 = YI - Y(JO)
        RX2 = XI - X(JP)
        RY2 = YI - Y(JP)
C
        SX = (X(JP) - X(JO)) * DSIO
        SY = (Y(JP) - Y(JO)) * DSIO
C
        X1 = SX*RX1 + SY*RY1
        X2 = SX*RX2 + SY*RY2
        YY = SX*RY1 - SY*RX1
C
        RS1 = RX1*RX1 + RY1*RY1
        RS2 = RX2*RX2 + RY2*RY2
C
        IF(IO.GE.N+1 .AND. IO.LE.N+NW) THEN
         SGN = 1.0
        ELSE
         SGN = SIGN(1.0,YY)
        ENDIF
C
        IF(IO.NE.JO .AND. RS1.GT.0.0) THEN
         G1 = LOG(RS1)
         T1 = ATAN2(SGN*X1,SGN*YY) - (0.5 - 0.5*SGN)*PI
        ELSE
         G1 = 0.0
         T1 = 0.0
        ENDIF
C
        IF(IO.NE.JP .AND. RS2.GT.0.0) THEN
         G2 = LOG(RS2)
         T2 = ATAN2(SGN*X2,SGN*YY) - (0.5 - 0.5*SGN)*PI
        ELSE
         G2 = 0.0
         T2 = 0.0
        ENDIF
C
        X1I = SX*NXI + SY*NYI
        X2I = SX*NXI + SY*NYI
        YYI = SX*NYI - SY*NXI
C
C------- set up midpoint quantities
         X0 = 0.5*(X1+X2)
         RS0 = X0*X0 + YY*YY
         G0 = LOG(RS0)
         T0 = ATAN2(SGN*X0,SGN*YY) - (0.5 - 0.5*SGN)*PI
C
C------- calculate source contribution to Psi  for  1-0  half-panel
         DXINV = 1.0/(X1-X0)
         PSUM = X0*(T0-APAN) - X1*(T1-APAN) + 0.5*YY*(G1-G0)
         PDIF = ((X1+X0)*PSUM + RS1*(T1-APAN) - RS0*(T0-APAN)
     &        + (X0-X1)*YY) * DXINV
C
         PSX1 =  -(T1-APAN)
         PSX0 =    T0-APAN
         PSYY =  0.5*(G1-G0)
C
         PDX1 = ((X1+X0)*PSX1 + PSUM + 2.0*X1*(T1-APAN) - PDIF) * DXINV
         PDX0 = ((X1+X0)*PSX0 + PSUM - 2.0*X0*(T0-APAN) + PDIF) * DXINV
         PDYY = ((X1+X0)*PSYY + 2.0*(X0-X1 + YY*(T1-T0))      ) * DXINV
C
         DSM = SQRT((X(JP)-X(JM))**2 + (Y(JP)-Y(JM))**2)
         DSIM = 1.0/DSM
C
CCC         SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC         SIG1 = (SIG(JP) - SIG(JM))*DSIM
CCC         SSUM = SIG0 + SIG1
CCC         SDIF = SIG0 - SIG1
C
         SSUM = (SIG(JP) - SIG(JO))*DSIO + (SIG(JP) - SIG(JM))*DSIM
         SDIF = (SIG(JP) - SIG(JO))*DSIO - (SIG(JP) - SIG(JM))*DSIM
C
         PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         DZDM(JM) = DZDM(JM) + QOPI*(-PSUM*DSIM + PDIF*DSIM)
         DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*DSIO - PDIF*DSIO)
         DZDM(JP) = DZDM(JP) + QOPI*( PSUM*(DSIO+DSIM)
     &                                          + PDIF*(DSIO-DSIM))
C
C------- dPsi/dni
         PSNI = PSX1*X1I + PSX0*(X1I+X2I)*0.5 + PSYY*YYI
         PDNI = PDX1*X1I + PDX0*(X1I+X2I)*0.5 + PDYY*YYI
         PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         DQDM(JM) = DQDM(JM) + QOPI*(-PSNI*DSIM + PDNI*DSIM)
         DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*DSIO - PDNI*DSIO)
         DQDM(JP) = DQDM(JP) + QOPI*( PSNI*(DSIO+DSIM)
     &                                          + PDNI*(DSIO-DSIM))
C
C
C------- calculate source contribution to Psi  for  0-2  half-panel
         DXINV = 1.0/(X0-X2)
         PSUM = X2*(T2-APAN) - X0*(T0-APAN) + 0.5*YY*(G0-G2)
         PDIF = ((X0+X2)*PSUM + RS0*(T0-APAN) - RS2*(T2-APAN)
     &        + (X2-X0)*YY) * DXINV
C
         PSX0 =  -(T0-APAN)
         PSX2 =    T2-APAN
         PSYY =  0.5*(G0-G2)
C
         PDX0 = ((X0+X2)*PSX0 + PSUM + 2.0*X0*(T0-APAN) - PDIF) * DXINV
         PDX2 = ((X0+X2)*PSX2 + PSUM - 2.0*X2*(T2-APAN) + PDIF) * DXINV
         PDYY = ((X0+X2)*PSYY + 2.0*(X2-X0 + YY*(T0-T2))      ) * DXINV
C
         DSP = SQRT((X(JQ)-X(JO))**2 + (Y(JQ)-Y(JO))**2)
         DSIP = 1.0/DSP
C
CCC         SIG2 = (SIG(JQ) - SIG(JO))*DSIP
CCC         SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC         SSUM = SIG2 + SIG0
CCC         SDIF = SIG2 - SIG0
C
         SSUM = (SIG(JQ) - SIG(JO))*DSIP + (SIG(JP) - SIG(JO))*DSIO
         SDIF = (SIG(JQ) - SIG(JO))*DSIP - (SIG(JP) - SIG(JO))*DSIO
C
         PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF)
C
C------- dPsi/dm
         DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*(DSIP+DSIO)
     &                                          - PDIF*(DSIP-DSIO))
         DZDM(JP) = DZDM(JP) + QOPI*( PSUM*DSIO - PDIF*DSIO)
         DZDM(JQ) = DZDM(JQ) + QOPI*( PSUM*DSIP + PDIF*DSIP)
C
C------- dPsi/dni
         PSNI = PSX0*(X1I+X2I)*0.5 + PSX2*X2I + PSYY*YYI
         PDNI = PDX0*(X1I+X2I)*0.5 + PDX2*X2I + PDYY*YYI
         PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*(DSIP+DSIO)
     &                                          - PDNI*(DSIP-DSIO))
         DQDM(JP) = DQDM(JP) + QOPI*( PSNI*DSIO - PDNI*DSIO)
         DQDM(JQ) = DQDM(JQ) + QOPI*( PSNI*DSIP + PDNI*DSIP)
C
   20 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Calculates source panel influence coefficient
C     matrix for current airfoil and wake geometry.
C
C===================================================================70
      SUBROUTINE QDCALC

      use my_equivalence, only : my_equiv_3_2
      use xfoil_inc
C
C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE)
     &  WRITE(*,*) 'Calculating source influence matrix ...'
C
      IF(.NOT.LADIJ) THEN
C
C----- calculate source influence matrix for airfoil surface if it doesn't exist
       DO 10 J=1, N
C
C------- multiply each dPsi/Sig vector by inverse of factored dPsi/dGam matrix
         CALL BAKSUB(IQX,N+1,AIJ,AIJPIV,BIJ(1,J))

C        DP mod: copy to VM; used to replace equivalence statement
         call my_equiv_3_2(VM, BIJ, (/ 1, 1, 1 /), (/ 1, 1 /),
     ?                    (/ 1, J /), 1)

C
C------- store resulting dGam/dSig = dQtan/dSig vector
         DO 105 I=1, N
           DIJ(I,J) = BIJ(I,J)
  105    CONTINUE
C
   10  CONTINUE
       LADIJ = .TRUE.
C
      ENDIF
C
C---- set up coefficient matrix of dPsi/dm on airfoil surface
      DO 20 I=1, N
        CALL PSWLIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_N)
        DO 202 J=N+1, N+NW
          BIJ(I,J) = -DZDM(J)

C         DP mod: copy to VM; used to replace equivalence statement
          call my_equiv_3_2(VM, BIJ, (/ 1, 1, 1 /), (/ 1, 1 /),
     ?                     (/ I, J /), 1)

  202   CONTINUE
   20 CONTINUE
C
C---- set up Kutta condition (no direct source influence)
      DO 32 J=N+1, N+NW
        BIJ(N+1,J) = 0.

C       DP mod: copy to VM; used to replace equivalence statement
        call my_equiv_3_2(VM, BIJ, (/ 1, 1, 1 /), (/ 1, 1 /),
     ?                   (/ N+1, J /), 1)

   32 CONTINUE
C
C---- sharp TE gamma extrapolation also has no source influence
      IF(SHARP) THEN
       DO 34 J=N+1, N+NW
         BIJ(N,J) = 0.

C        DP mod: copy to VM; used to replace equivalence statement
         call my_equiv_3_2(VM, BIJ, (/ 1, 1, 1 /), (/ 1, 1 /),
     ?                    (/ N, J /), 1)

   34  CONTINUE
      ENDIF
C
C---- multiply by inverse of factored dPsi/dGam matrix
      DO 40 J=N+1, N+NW
        CALL BAKSUB(IQX,N+1,AIJ,AIJPIV,BIJ(1,J))

C       DP mod: copy to VM; used to replace equivalence statement
        call my_equiv_3_2(VM, BIJ, (/ 1, 1, 1 /), (/ 1, 1 /),
     ?                   (/ 1, J /), 1)

   40 CONTINUE
C
C---- set the source influence matrix for the wake sources
      DO 50 I=1, N
        DO 510 J=N+1, N+NW
          DIJ(I,J) = BIJ(I,J)
  510   CONTINUE
   50 CONTINUE
C
C**** Now we need to calculate the influence of sources on the wake velocities
C
C---- calculcate dQtan/dGam and dQtan/dSig at the wake points
      DO 70 I=N+1, N+NW
C
        IW = I-N
C
C------ airfoil contribution at wake panel node
        CALL PSILIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_N,.FALSE.,.TRUE.)
C
        DO 710 J=1, N
          CIJ(IW,J) = DQDG(J)
C
C         DP mod: copy to VM; used to replace equivalence statement
          call my_equiv_3_2(VM, CIJ, (/ 1, 1, IZX/2 /), (/ 1, 1 /),
     ?                     (/ IW, J /), 1)
C
  710   CONTINUE
C  
        DO 720 J=1, N
          DIJ(I,J) = DQDM(J)
  720   CONTINUE
C
C------ wake contribution
        CALL PSWLIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_N)
C
        DO 730 J=N+1, N+NW
          DIJ(I,J) = DQDM(J)
  730   CONTINUE
C
   70 CONTINUE
C
C---- add on effect of all sources on airfoil vorticity which effects wake Qtan
      DO 80 I=N+1, N+NW
        IW = I-N
C
C------ airfoil surface source contribution first
        DO 810 J=1, N
          SUM = 0.
          DO 8100 K=1, N
            SUM = SUM + CIJ(IW,K)*DIJ(K,J)
 8100     CONTINUE
          DIJ(I,J) = DIJ(I,J) + SUM
  810   CONTINUE
C
C------ wake source contribution next
        DO 820 J=N+1, N+NW
          SUM = 0.
          DO 8200 K=1, N
            SUM = SUM + CIJ(IW,K)*BIJ(K,J)
 8200     CONTINUE
          DIJ(I,J) = DIJ(I,J) + SUM
  820   CONTINUE
C
   80 CONTINUE
C
C---- make sure first wake point has same velocity as trailing edge
      DO 90 J=1, N+NW
        DIJ(N+1,J) = DIJ(N,J)
   90 CONTINUE
C
      LWDIJ = .TRUE.
C
      RETURN
      END

C===================================================================70
C
C     Calculates the "inverse" spline function S(X).    |
C     Since S(X) can be multi-valued or not defined,    |
C     this is not a "black-box" routine.  The calling   |
C     program must pass via SI a sufficiently good      |
C     initial guess for S(XI).                          |
C                                                       |
C     XI      specified X value       (input)           |
C     SI      calculated S(XI) value  (input,output)    |
C     X,XS,S  usual spline arrays     (input)           |
C                                                       |
C
C===================================================================70
      SUBROUTINE SINVRT(SI,XI,X,XS,S,N,SILENT_MODE)

      DIMENSION X(N), XS(N), S(N)
      LOGICAL :: SILENT_MODE
C
      SISAV = SI
C
      DO 10 ITER=1, 10
        RES  = SEVAL(SI,X,XS,S,N) - XI
        RESP = DEVAL(SI,X,XS,S,N)
        DS = -RES/RESP
        SI = SI + DS
        IF(ABS(DS/(S(N)-S(1))) .LT. 1.0E-5) RETURN
   10 CONTINUE
C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE) THEN
        WRITE(*,*)
     &    'SINVRT: spline inversion failed. Input value returned.'
      ENDIF
      SI = SISAV
C
      RETURN
      END ! SINVRT

C===================================================================70
C
C Computes top and bottom hinge locations for flap
C
C===================================================================70
      SUBROUTINE GETXYF(X,XP,Y,YP,S,N, TOPS,BOTS,XF,YF,SILENT_MODE)
      DIMENSION X(N),XP(N),Y(N),YP(N),S(N)
C
C      DP mod: XF is specified
C      IF(XF .EQ. -999.0)
C     &  CALL ASKR('Enter flap hinge x location^',XF)
C
C---- find top and bottom y at hinge x location
      TOPS = S(1) + (X(1) - XF)
      BOTS = S(N) - (X(N) - XF)
C     DP mod: added SILENT_MODE option
      CALL SINVRT(TOPS,XF,X,XP,S,N,SILENT_MODE)      
      CALL SINVRT(BOTS,XF,X,XP,S,N,SILENT_MODE)      
      TOPY = SEVAL(TOPS,Y,YP,S,N)
      BOTY = SEVAL(BOTS,Y,YP,S,N)
C
C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE) THEN
        WRITE(*,1000) TOPY, BOTY
      ENDIF
 1000 FORMAT(/'  Top    surface:  y =', F8.4,'     y/t = 1.0'
     &       /'  Bottom surface:  y =', F8.4,'     y/t = 0.0')
C
C      DP mod: YF is specified
C      IF(YF .EQ. -999.0)
C     & CALL ASKR(
C     &  'Enter flap hinge y location (or 999 to specify y/t)^',YF)
CC
C      IF(YF .EQ. 999.0) THEN
C        CALL ASKR('Enter flap hinge relative y/t location^',YREL)
C        YF = TOPY*YREL + BOTY*(1.0-YREL)
C      ENDIF
C
      RETURN
      END ! GETXYF

C===================================================================70
C
C     Returns .TRUE. if point XF,YF 
C     is inside contour X(i),Y(i).
C
C===================================================================70
      LOGICAL FUNCTION INSIDE(X,Y,N, XF,YF)
      DIMENSION X(N),Y(N)
C
C---- integrate subtended angle around airfoil perimeter
      ANGLE = 0.0
      DO 10 I=1, N
        IP = I+1
        IF(I.EQ.N) IP = 1
        XB1 = X(I)  - XF
        YB1 = Y(I)  - YF
        XB2 = X(IP) - XF
        YB2 = Y(IP) - YF
        ANGLE = ANGLE + (XB1*YB2 - YB1*XB2)
     &                   / SQRT((XB1**2 + YB1**2)*(XB2**2 + YB2**2))
 10   CONTINUE
C
C---- angle = 0 if XF,YF is outside, angle = +/- 2 pi  if XF,YF is inside
      INSIDE = ABS(ANGLE) .GT. 1.0
C
      RETURN
      END ! INSIDE

C===================================================================70
C
C     Returns arc length points S1,S2 at flap surface break
C     locations.  S1 is on fixed airfoil part, S2 is on flap.
C     The points are defined according to two cases:
C
C
C     If DEL > 0:  Surface will be eliminated in S1 < s < S2
C
C     Returns the arc length values S1,S2 of the endpoints
C     of the airfoil surface segment which "disappears" as a
C     result of the flap deflection.  The line segments between
C     these enpoints and the flap hinge point (XBF,YBF) have
C     an included angle of DEL.  DEL is therefore the flap
C     deflection which will join up the points at S1,S2.
C     SS is an approximate arc length value near S1 and S2.
C     It is used as an initial guess for the Newton loop 
C     for S1 and S2.
C
C
C     If DEL = 0:  Surface will be created at s = S1 = S2
C
C     If DEL=0, then S1,S2 will cooincide, and will be located
C     on the airfoil surface where the segment joining the
C     point at S1,S2 and the hinge point is perpendicular to
C     the airfoil surface.  This will be the point where the
C     airfoil surface must be broken to permit a gap to open
C     as a result of the flap deflection.
C
C===================================================================70
      SUBROUTINE SSS(SS,S1,S2,DEL,XBF,YBF,X,XP,Y,YP,S,N,ISIDE,
     &               SILENT_MODE)
      DIMENSION X(*),XP(*),Y(*),YP(*),S(*)
C
C---- convergence epsilon
      DATA EPS / 1.0E-5 /
C
      STOT = ABS( S(N) - S(1) )
C
      SIND = SIN(0.5*ABS(DEL))
C
      SSGN = 1.0
      IF(ISIDE.EQ.1) SSGN = -1.0
C
C---- initial guesses for S1, S2
      RSQ = (SEVAL(SS,X,XP,S,N)-XBF)**2 + (SEVAL(SS,Y,YP,S,N)-YBF)**2
      S1 = SS - (SIND*SQRT(RSQ) + EPS*STOT)*SSGN
      S2 = SS + (SIND*SQRT(RSQ) + EPS*STOT)*SSGN
C
C---- Newton iteration loop
      DO 10 ITER=1, 10
        X1  = SEVAL(S1,X,XP,S,N)
        X1P = DEVAL(S1,X,XP,S,N)
        Y1  = SEVAL(S1,Y,YP,S,N)
        Y1P = DEVAL(S1,Y,YP,S,N)
C
        X2  = SEVAL(S2,X,XP,S,N)
        X2P = DEVAL(S2,X,XP,S,N)
        Y2  = SEVAL(S2,Y,YP,S,N)
        Y2P = DEVAL(S2,Y,YP,S,N)
C
        R1SQ = (X1-XBF)**2 + (Y1-YBF)**2
        R2SQ = (X2-XBF)**2 + (Y2-YBF)**2
        R1 = SQRT(R1SQ)
        R2 = SQRT(R2SQ)
C
        RRSQ = (X1-X2)**2 + (Y1-Y2)**2
        RR = SQRT(RRSQ)
C
        IF(R1.LE.EPS*STOT .OR. R2.LE.EPS*STOT) THEN
         S1 = SS
         S2 = SS
         RETURN
        ENDIF
C
        R1_S1 = (X1P*(X1-XBF) + Y1P*(Y1-YBF))/R1
        R2_S2 = (X2P*(X2-XBF) + Y2P*(Y2-YBF))/R2
C
        IF(SIND.GT.0.01) THEN
C
         IF(RR.EQ.0.0) RETURN
C
         RR_S1 =  (X1P*(X1-X2) + Y1P*(Y1-Y2))/RR
         RR_S2 = -(X2P*(X1-X2) + Y2P*(Y1-Y2))/RR
C
C------- Residual 1: set included angle via dot product
         RS1 = ((XBF-X1)*(X2-X1) + (YBF-Y1)*(Y2-Y1))/RR - SIND*R1
         A11 = ((XBF-X1)*( -X1P) + (YBF-Y1)*( -Y1P))/RR
     &       + ((  -X1P)*(X2-X1) + (  -Y1P)*(Y2-Y1))/RR
     &       - ((XBF-X1)*(X2-X1) + (YBF-Y1)*(Y2-Y1))*RR_S1/RRSQ
     &       - SIND*R1_S1
         A12 = ((XBF-X1)*(X2P  ) + (YBF-Y1)*(Y2P  ))/RR
     &       - ((XBF-X1)*(X2-X1) + (YBF-Y1)*(Y2-Y1))*RR_S2/RRSQ
C
C------- Residual 2: set equal length segments
         RS2 = R1 - R2
         A21 = R1_S1
         A22 =    - R2_S2
C
        ELSE
C
C------- Residual 1: set included angle via small angle approximation
         RS1 = (R1+R2)*SIND + (S1 - S2)*SSGN
         A11 =  R1_S1 *SIND + SSGN
         A12 =  R2_S2 *SIND - SSGN
C
C------- Residual 2: set vector sum of line segments beteen the 
C-       endpoints and flap hinge to be perpendicular to airfoil surface.
         X1PP = D2VAL(S1,X,XP,S,N)
         Y1PP = D2VAL(S1,Y,YP,S,N)
         X2PP = D2VAL(S2,X,XP,S,N)
         Y2PP = D2VAL(S2,Y,YP,S,N)
C
         XTOT = X1+X2 - 2.0*XBF
         YTOT = Y1+Y2 - 2.0*YBF
C
         RS2 = XTOT*(X1P+X2P) + YTOT*(Y1P+Y2P)
         A21 =  X1P*(X1P+X2P) +  Y1P*(Y1P+Y2P) + XTOT*X1PP + YTOT*Y1PP
         A22 =  X2P*(X1P+X2P) +  Y2P*(Y1P+Y2P) + XTOT*X2PP + YTOT*Y2PP
C
        ENDIF
C
        DET =   A11*A22 - A12*A21
        DS1 = -(RS1*A22 - A12*RS2) / DET
        DS2 = -(A11*RS2 - RS1*A21) / DET
C
        DS1 = MIN( DS1 , 0.01*STOT )
        DS1 = MAX( DS1 , -.01*STOT )
        DS2 = MIN( DS2 , 0.01*STOT )
        DS2 = MAX( DS2 , -.01*STOT )
C
        S1 = S1 + DS1
        S2 = S2 + DS2
        IF(ABS(DS1)+ABS(DS2) .LT. EPS*STOT ) GO TO 11
   10 CONTINUE
C     DP mod: added option for SILENT_MODE
      IF (.NOT. SILENT_MODE) THEN
        WRITE(*,*) 'SSS: failed to converge subtending angle points'
      ENDIF
      S1 = SS
      S2 = SS
C
   11 CONTINUE
C
C---- make sure points are identical if included angle is zero.
      IF(DEL.EQ.0.0) THEN
       S1 = 0.5*(S1+S2)
       S2 = S1
      ENDIF
C
      RETURN
      END

C===================================================================70
C
C     Modifies buffer airfoil for a deflected flap.
C     Points may be added/subtracted in the flap
C     break vicinity to clean things up.
C
C===================================================================70
      SUBROUTINE FLAP(XBF,YBF,DDEF,SILENT_MODE)

C FIXME: Before using, make sure XB, XBP, YB, YBP, SB, NB are set
C properly from airfoil generated by PANGEN.
      use xfoil_inc

      LOGICAL LCHANGE
      DIMENSION RINPUT(*)
C
      LOGICAL INSID
      LOGICAL INSIDE
      LOGICAL LT1NEW,LT2NEW,LB1NEW,LB2NEW
C
C      DP mod: below is used for plotting
C      SHT = CH * MAX(XSF,YSF)
C
C      DP mod: XBF and YBF changed to separate inputs
C      IF(NINPUT.GE.2) THEN
C       XBF = RINPUT(1)
C       YBF = RINPUT(2)
C      ELSE
C       XBF = -999.0
C       YBF = -999.0
C      ENDIF
C
C     DP mod: added SILENT_MODE option
      CALL GETXYF(XB,XBP,YB,YBP,SB,NB, TOPS,BOTS,XBF,YBF,SILENT_MODE)
      INSID = INSIDE(XB,YB,NB,XBF,YBF)
C
C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE) THEN
        WRITE(*,1050) XBF, YBF
      ENDIF
 1050 FORMAT(/' Flap hinge: x,y =', 2F9.5 )
C
C     DP mod: DDEF is specified
C      IF(NINPUT.GE.3) THEN
C       DDEF = RINPUT(3)
C      ELSE
C       DDEF = 0.0
C       CALL ASKR('Enter flap deflection in degrees (+ down)^',DDEF)
C      ENDIF
      RDEF = DDEF*PI/180.0
      IF(RDEF .EQ. 0.0) RETURN
C
C
      IF(INSID) THEN
        ATOP = MAX( 0.0 , -RDEF )
        ABOT = MAX( 0.0 ,  RDEF )
      ELSE
        CHX = DEVAL(BOTS,XB,XBP,SB,NB) - DEVAL(TOPS,XB,XBP,SB,NB)
        CHY = DEVAL(BOTS,YB,YBP,SB,NB) - DEVAL(TOPS,YB,YBP,SB,NB)
        FVX = SEVAL(BOTS,XB,XBP,SB,NB) + SEVAL(TOPS,XB,XBP,SB,NB)
        FVY = SEVAL(BOTS,YB,YBP,SB,NB) + SEVAL(TOPS,YB,YBP,SB,NB)
        CRSP = CHX*(YBF-0.5*FVY) - CHY*(XBF-0.5*FVX)
        IF(CRSP .GT. 0.0) THEN
C-------- flap hinge is above airfoil
          ATOP = MAX( 0.0 ,  RDEF )
          ABOT = MAX( 0.0 ,  RDEF )
        ELSE
C-------- flap hinge is below airfoil
          ATOP = MAX( 0.0 , -RDEF )
          ABOT = MAX( 0.0 , -RDEF )
        ENDIF
      ENDIF
C
C---- find upper and lower surface break arc length values...
C     DP mod: added option for SILENT_MODE
      CALL SSS(TOPS,ST1,ST2,ATOP,XBF,YBF,XB,XBP,YB,YBP,SB,NB,1,
     &         SILENT_MODE)
      CALL SSS(BOTS,SB1,SB2,ABOT,XBF,YBF,XB,XBP,YB,YBP,SB,NB,2,
     &         SILENT_MODE)
C
C---- ... and x,y coordinates
      XT1 = SEVAL(ST1,XB,XBP,SB,NB)
      YT1 = SEVAL(ST1,YB,YBP,SB,NB)
      XT2 = SEVAL(ST2,XB,XBP,SB,NB)
      YT2 = SEVAL(ST2,YB,YBP,SB,NB)
      XB1 = SEVAL(SB1,XB,XBP,SB,NB)
      YB1 = SEVAL(SB1,YB,YBP,SB,NB)
      XB2 = SEVAL(SB2,XB,XBP,SB,NB)
      YB2 = SEVAL(SB2,YB,YBP,SB,NB)
C
C
C     DP mod: added option for SILENT_MODE
      IF (.NOT. SILENT_MODE) THEN
        WRITE(*,1100) XT1, YT1, XT2, YT2,
     &                XB1, YB1, XB2, YB2
      ENDIF
 1100 FORMAT(/' Top breaks: x,y =  ', 2F9.5, 4X, 2F9.5
     &       /' Bot breaks: x,y =  ', 2F9.5, 4X, 2F9.5)
C
C---- find points adjacent to breaks
      DO 5 I=1, NB-1
        IF(SB(I).LE.ST1 .AND. SB(I+1).GT.ST1) IT1 = I+1
        IF(SB(I).LT.ST2 .AND. SB(I+1).GE.ST2) IT2 = I
        IF(SB(I).LE.SB1 .AND. SB(I+1).GT.SB1) IB1 = I
        IF(SB(I).LT.SB2 .AND. SB(I+1).GE.SB2) IB2 = I+1
    5 CONTINUE
C
      DSAVG = (SB(NB)-SB(1))/FLOAT(NB-1)
C
C---- smallest fraction of s increments i+1 and i+2 away from break point
      SFRAC = 0.33333
C
      IF(ATOP .NE. 0.0) THEN
        ST1P = ST1 + SFRAC*(SB(IT1  )-ST1)
        ST1Q = ST1 + SFRAC*(SB(IT1+1)-ST1)
        IF(SB(IT1) .LT. ST1Q) THEN
C-------- simply move adjacent point to ideal SFRAC location
          XT1NEW = SEVAL(ST1Q,XB,XBP,SB,NB)
          YT1NEW = SEVAL(ST1Q,YB,YBP,SB,NB)
          LT1NEW = .FALSE.
        ELSE
C-------- make new point at SFRAC location
          XT1NEW = SEVAL(ST1P,XB,XBP,SB,NB)
          YT1NEW = SEVAL(ST1P,YB,YBP,SB,NB)
          LT1NEW = .TRUE.
        ENDIF
C
        ST2P = ST2 + SFRAC*(SB(IT2 )-ST2)
        IT2Q = MAX(IT2-1,1)
        ST2Q = ST2 + SFRAC*(SB(IT2Q)-ST2)
        IF(SB(IT2) .GT. ST2Q) THEN
C-------- simply move adjacent point
          XT2NEW = SEVAL(ST2Q,XB,XBP,SB,NB)
          YT2NEW = SEVAL(ST2Q,YB,YBP,SB,NB)
          LT2NEW = .FALSE.
        ELSE
C-------- make new point
          XT2NEW = SEVAL(ST2P,XB,XBP,SB,NB)
          YT2NEW = SEVAL(ST2P,YB,YBP,SB,NB)
          LT2NEW = .TRUE.
        ENDIF
      ENDIF
C
      IF(ABOT .NE. 0.0) THEN
        SB1P = SB1 + SFRAC*(SB(IB1  )-SB1)
        SB1Q = SB1 + SFRAC*(SB(IB1-1)-SB1)
        IF(SB(IB1) .GT. SB1Q) THEN
C-------- simply move adjacent point
          XB1NEW = SEVAL(SB1Q,XB,XBP,SB,NB)
          YB1NEW = SEVAL(SB1Q,YB,YBP,SB,NB)
          LB1NEW = .FALSE.
        ELSE
C-------- make new point
          XB1NEW = SEVAL(SB1P,XB,XBP,SB,NB)
          YB1NEW = SEVAL(SB1P,YB,YBP,SB,NB)
          LB1NEW = .TRUE.
        ENDIF
C
        SB2P = SB2 + SFRAC*(SB(IB2 )-SB2)
        IB2Q = MIN(IB2+1,NB)
        SB2Q = SB2 + SFRAC*(SB(IB2Q)-SB2)
        IF(SB(IB2) .LT. SB2Q) THEN
C-------- simply move adjacent point
          XB2NEW = SEVAL(SB2Q,XB,XBP,SB,NB)
          YB2NEW = SEVAL(SB2Q,YB,YBP,SB,NB)
          LB2NEW = .FALSE.
        ELSE
C-------- make new point
          XB2NEW = SEVAL(SB2P,XB,XBP,SB,NB)
          YB2NEW = SEVAL(SB2P,YB,YBP,SB,NB)
          LB2NEW = .TRUE.
        ENDIF
      ENDIF
C
cc      DSTOP = ABS(SB(IT2)-SB(IT1))
cc      DSBOT = ABS(SB(IB2)-SB(IB1))
C
      SIND = SIN(RDEF)
      COSD = COS(RDEF)
C
C---- rotate flap points about the hinge point (XBF,YBF)
      DO 10 I=1, NB
        IF(I.GE.IT1 .AND. I.LE.IB1) GO TO 10
C
        XBAR = XB(I) - XBF
        YBAR = YB(I) - YBF
C
        XB(I) = XBF  +  XBAR*COSD  +  YBAR*SIND
        YB(I) = YBF  -  XBAR*SIND  +  YBAR*COSD
   10 CONTINUE
C
      IDIF = IT1-IT2-1
      IF(IDIF.GT.0) THEN
C----- delete points on upper airfoil surface which "disappeared".
       NB  = NB -IDIF
       IT1 = IT1-IDIF
       IB1 = IB1-IDIF
       IB2 = IB2-IDIF
       DO 21 I=IT2+1, NB
         SB(I) = SB(I+IDIF)
         XB(I) = XB(I+IDIF)
         YB(I) = YB(I+IDIF)
   21  CONTINUE
      ENDIF
C
      IDIF = IB2-IB1-1
      IF(IDIF.GT.0) THEN
C----- delete points on lower airfoil surface which "disappeared".
       NB  = NB -IDIF
       IB2 = IB2-IDIF
       DO 22 I=IB1+1, NB
         SB(I) = SB(I+IDIF)
         XB(I) = XB(I+IDIF)
         YB(I) = YB(I+IDIF)
   22  CONTINUE
      ENDIF
C
C
      IF(ATOP .EQ. 0.0) THEN
C
C------ arc length of newly created surface on top of airfoil
        DSNEW = ABS(RDEF)*SQRT((XT1-XBF)**2 + (YT1-YBF)**2)
C
C------ number of points to be added to define newly created surface
        NPADD = INT(1.5*DSNEW/DSAVG + 1.0)
ccc     NPADD = INT(1.5*DSNEW/DSTOP + 1.0)
C
C------ skip everything if no points are to be added
        IF(NPADD.EQ.0) GO TO 35
C
C------ increase coordinate array length to make room for the new point(s)
        NB  = NB +NPADD
        IT1 = IT1+NPADD
        IB1 = IB1+NPADD
        IB2 = IB2+NPADD
        DO 30 I=NB, IT1, -1
          XB(I) = XB(I-NPADD)
          YB(I) = YB(I-NPADD)
   30   CONTINUE
C
C------ add new points along the new surface circular arc segment
        DANG = RDEF / FLOAT(NPADD)
        XBAR = XT1 - XBF
        YBAR = YT1 - YBF
        DO 31 IP=1, NPADD
          ANG = DANG*(FLOAT(IP) - 0.5)
          CA = COS(ANG)
          SA = SIN(ANG)
C
          XB(IT1-IP) = XBF  +  XBAR*CA + YBAR*SA
          YB(IT1-IP) = YBF  -  XBAR*SA + YBAR*CA
   31   CONTINUE
C
      ELSE
C
C------ set point in the corner and possibly two adjacent points
        NPADD = 1
        IF(LT2NEW) NPADD = NPADD+1
        IF(LT1NEW) NPADD = NPADD+1
C
        NB  = NB +NPADD
        IT1 = IT1+NPADD
        IB1 = IB1+NPADD
        IB2 = IB2+NPADD
        DO 33 I=NB, IT1, -1
          XB(I) = XB(I-NPADD)
          YB(I) = YB(I-NPADD)
   33   CONTINUE
C
        IF(LT1NEW) THEN
         XB(IT1-1) = XT1NEW
         YB(IT1-1) = YT1NEW
         XB(IT1-2) = XT1
         YB(IT1-2) = YT1
        ELSE
         XB(IT1  ) = XT1NEW
         YB(IT1  ) = YT1NEW
         XB(IT1-1) = XT1
         YB(IT1-1) = YT1
        ENDIF
C
        XBAR = XT2NEW - XBF
        YBAR = YT2NEW - YBF
        IF(LT2NEW) THEN
          XB(IT2+1) = XBF  +  XBAR*COSD + YBAR*SIND
          YB(IT2+1) = YBF  -  XBAR*SIND + YBAR*COSD
        ELSE
          XB(IT2  ) = XBF  +  XBAR*COSD + YBAR*SIND
          YB(IT2  ) = YBF  -  XBAR*SIND + YBAR*COSD
        ENDIF
C
      ENDIF
   35 CONTINUE
C
C
      IF(ABOT .EQ. 0.0) THEN
C
C------ arc length of newly created surface on top of airfoil
        DSNEW = ABS(RDEF)*SQRT((XB1-XBF)**2 + (YB1-YBF)**2)
C
C------ number of points to be added to define newly created surface
        NPADD = INT(1.5*DSNEW/DSAVG + 1.0)
ccc     NPADD = INT(1.5*DSNEW/DSBOT + 1.0)
C
C------ skip everything if no points are to be added
        IF(NPADD.EQ.0) GO TO 45
C
C------ increase coordinate array length to make room for the new point(s)
        NB  = NB +NPADD
        IB2 = IB2+NPADD
        DO 40 I=NB, IB2, -1
          XB(I) = XB(I-NPADD)
          YB(I) = YB(I-NPADD)
   40   CONTINUE
C
C------ add new points along the new surface circular arc segment
        DANG = RDEF / FLOAT(NPADD)
        XBAR = XB1 - XBF
        YBAR = YB1 - YBF
        DO 41 IP=1, NPADD
          ANG = DANG*(FLOAT(IP) - 0.5)
          CA = COS(ANG)
          SA = SIN(ANG)
C
          XB(IB1+IP) = XBF  +  XBAR*CA + YBAR*SA
          YB(IB1+IP) = YBF  -  XBAR*SA + YBAR*CA
   41   CONTINUE
C
      ELSE

C------ set point in the corner and possibly two adjacent points
        NPADD = 1
        IF(LB2NEW) NPADD = NPADD+1
        IF(LB1NEW) NPADD = NPADD+1
C
        NB  = NB +NPADD
        IB2 = IB2+NPADD
        DO 43 I=NB, IB2, -1
          XB(I) = XB(I-NPADD)
          YB(I) = YB(I-NPADD)
   43   CONTINUE
C
        IF(LB1NEW) THEN
         XB(IB1+1) = XB1NEW
         YB(IB1+1) = YB1NEW
         XB(IB1+2) = XB1
         YB(IB1+2) = YB1
        ELSE
         XB(IB1  ) = XB1NEW
         YB(IB1  ) = YB1NEW
         XB(IB1+1) = XB1
         YB(IB1+1) = YB1
        ENDIF
C
        XBAR = XB2NEW - XBF
        YBAR = YB2NEW - YBF
        IF(LB2NEW) THEN
          XB(IB2-1) = XBF  +  XBAR*COSD + YBAR*SIND
          YB(IB2-1) = YBF  -  XBAR*SIND + YBAR*COSD
        ELSE
          XB(IB2  ) = XBF  +  XBAR*COSD + YBAR*SIND
          YB(IB2  ) = YBF  -  XBAR*SIND + YBAR*COSD
        ENDIF
C
      ENDIF
   45 CONTINUE
C
      LGSAME = .FALSE.
C
C
C     DP note: got here
C---- check new geometry for splinter segments 
      STOL = 0.2
      CALL SCHECK(XB,YB,NB, STOL, LCHANGE)
C
C---- spline new geometry
      CALL SCALC(XB,YB,SB,NB)
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
C
      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
     &            EI11BA,EI22BA,APX1BA,APX2BA,
     &            EI11BT,EI22BT,APX1BT,APX2BT,
     &            THICKB,CAMBRB )
C
      LBFLAP = .TRUE.
C
      IF(LGSYM) THEN
       WRITE(*,*)
       WRITE(*,*) 'Disabling symmetry enforcement'
       LGSYM = .FALSE.
      ENDIF
C
C
      IF(.NOT.LPLOT) THEN
       CALL PLTINI
      ENDIF
C
C---- save current color and set new color
      CALL GETCOLOR(ICOL0)
C
      CALL NEWCOLORNAME('green')
      CALL PLOT((XBF-XOFF)*XSF,(YBF-YOFF)*YSF,3)
      CALL PLOT((XT1-XOFF)*XSF,(YT1-YOFF)*YSF,2)
      CALL PLOT((XBF-XOFF)*XSF,(YBF-YOFF)*YSF,3)
      CALL PLOT((XB1-XOFF)*XSF,(YB1-YOFF)*YSF,2)
C
      IF(ATOP .EQ. 0.0) THEN
        XBAR = XT1 - XBF
        YBAR = YT1 - YBF
        XT1C = XBF  +  XBAR*COSD + YBAR*SIND
        YT1C = YBF  -  XBAR*SIND + YBAR*COSD
        CALL PLOT((XBF -XOFF)*XSF,(YBF -YOFF)*YSF,3)
        CALL PLOT((XT1C-XOFF)*XSF,(YT1C-YOFF)*YSF,2)
      ENDIF
C
      IF(ABOT .EQ. 0.0) THEN
        XBAR = XB1 - XBF
        YBAR = YB1 - YBF
        XB1C = XBF  +  XBAR*COSD + YBAR*SIND
        YB1C = YBF  -  XBAR*SIND + YBAR*COSD
        CALL PLOT((XBF -XOFF)*XSF,(YBF -YOFF)*YSF,3)
        CALL PLOT((XB1C-XOFF)*XSF,(YB1C-YOFF)*YSF,2)
      ENDIF
C
      CALL NEWCOLORNAME('red')
      CALL PLSYMB((XBF-XOFF)*XSF,(YBF-YOFF)*YSF,0.5*SHT,1,0.0,0)
C
      CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
      CALL PLNEWP('magenta')
C
      LGEOPL = .FALSE.
C
      CALL NEWCOLOR(ICOL0)
      RETURN
      END ! FLAP
