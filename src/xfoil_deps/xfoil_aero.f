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

C  Copyright (C) 2014 -- 2016 Daniel Prosser (this modified version of 
C  XFoil code)
C  Original copyright (C) 2000 Mark Drela (original XFoil code)

C===================================================================70
C===================================================================70
      SUBROUTINE BLPINI

      use blpar_inc
C     
      SCCON = 5.6
      GACON = 6.70
      GBCON = 0.75
      GCCON = 18.0
      DLCON =  0.9
C
      CTRCON = 1.8
      CTRCEX = 3.3
C
      DUXCON = 1.0
C
      CTCON = 0.5/(GACON**2 * GBCON)
C
      CFFAC = 1.0
C
      RETURN
      END

C===================================================================70
C
C     Calculates current streamfunction Psi at panel node or wake node
C     I due to freestream and all bound vorticity Gam on the airfoil. 
C     Sensitivities of Psi with respect to alpha (Z_ALFA) and inverse
C     Qspec DOFs (Z_QDOF0,Z_QDOF1) which influence Gam in inverse cases.
C     Also calculates the sensitivity vector dPsi/dGam (DZDG).
C
C     If SIGLIN=True, then Psi includes the effects of the viscous
C     source distribution Sig and the sensitivity vector dPsi/dSig
C     (DZDM) is calculated.
C
C     If GEOLIN=True, then the geometric sensitivity vector dPsi/dn
C     is calculated, where n is the normal motion of the jth node.
C
C          Airfoil:  1   < I < N
C          Wake:     N+1 < I < N+NW
C
C===================================================================70
      SUBROUTINE PSILIN(I,XI,YI,NXI,NYI,PSI,PSI_NI,GEOLIN,SIGLIN)

      use xfoil_inc

      REAL*8 NXO, NYO, NXP, NYP, NXI, NYI
      LOGICAL GEOLIN,SIGLIN

C
C---- distance tolerance for determining if two points are the same
      SEPS = (S(N)-S(1)) * 1.0E-5
C
      IO = I
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
C
      DO 3 JO=1, N
        DZDG(JO) = 0.0
        DZDN(JO) = 0.0
        DQDG(JO) = 0.0
    3 CONTINUE
C
      DO 4 JO=1, N
        DZDM(JO) = 0.0
        DQDM(JO) = 0.0
    4 CONTINUE
C
      Z_QINF = 0.
      Z_ALFA = 0.
      Z_QDOF0 = 0.
      Z_QDOF1 = 0.
      Z_QDOF2 = 0.
      Z_QDOF3 = 0.
C
      PSI    = 0.
      PSI_NI = 0.
C
      QTAN1 = 0.
      QTAN2 = 0.
      QTANM = 0.
C
      IF(SHARP) THEN
       SCS = 1.0
       SDS = 0.0
      ELSE
       SCS = ANTE/DSTE
       SDS = ASTE/DSTE
      ENDIF
C
      DO 10 JO=1, N
        JP = JO+1
C
        JM = JO-1
        JQ = JP+1
C
        IF(JO.EQ.1) THEN
         JM = JO
        ELSE IF(JO.EQ.N-1) THEN
         JQ = JP
        ELSE IF(JO.EQ.N) THEN
         JP = 1
         IF((X(JO)-X(JP))**2 + (Y(JO)-Y(JP))**2 .LT. SEPS**2) GO TO 12
        ENDIF
C
        DSO = SQRT((X(JO)-X(JP))**2 + (Y(JO)-Y(JP))**2)
C
C------ skip null panel
        IF(DSO .EQ. 0.0) GO TO 10
C
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
C------ set reflection flag SGN to avoid branch problems with arctan
        IF(IO.GE.1 .AND. IO.LE.N) THEN
C------- no problem on airfoil surface
         SGN = 1.0
        ELSE
C------- make sure arctan falls between  -/+  Pi/2
         SGN = SIGN(1.0,YY)
        ENDIF
C
C------ set log(r^2) and arctan(x/y), correcting for reflection if any
        IF(IO.NE.JO .AND. RS1.GT.0.0) THEN
         G1 = LOG(RS1)
         T1 = ATAN2(SGN*X1,SGN*YY) + (0.5 - 0.5*SGN)*PI
        ELSE
         G1 = 0.0
         T1 = 0.0
        ENDIF
C
        IF(IO.NE.JP .AND. RS2.GT.0.0) THEN
         G2 = LOG(RS2)
         T2 = ATAN2(SGN*X2,SGN*YY) + (0.5 - 0.5*SGN)*PI
        ELSE
         G2 = 0.0
         T2 = 0.0
        ENDIF
C
        X1I = SX*NXI + SY*NYI
        X2I = SX*NXI + SY*NYI
        YYI = SX*NYI - SY*NXI
C
        IF(GEOLIN) THEN
         NXO = NX(JO)
         NYO = NY(JO)
         NXP = NX(JP)
         NYP = NY(JP)
C
         X1O =-((RX1-X1*SX)*NXO + (RY1-X1*SY)*NYO)*DSIO-(SX*NXO+SY*NYO)
         X1P = ((RX1-X1*SX)*NXP + (RY1-X1*SY)*NYP)*DSIO
         X2O =-((RX2-X2*SX)*NXO + (RY2-X2*SY)*NYO)*DSIO
         X2P = ((RX2-X2*SX)*NXP + (RY2-X2*SY)*NYP)*DSIO-(SX*NXP+SY*NYP)
         YYO = ((RX1+X1*SY)*NYO - (RY1-X1*SX)*NXO)*DSIO-(SX*NYO-SY*NXO)
         YYP =-((RX1-X1*SY)*NYP - (RY1+X1*SX)*NXP)*DSIO
        ENDIF
C
        IF(JO.EQ.N) GO TO 11
C
        IF(SIGLIN) THEN
C
C------- set up midpoint quantities
         X0 = 0.5*(X1+X2)
         RS0 = X0*X0 + YY*YY
         G0 = LOG(RS0)
         T0 = ATAN2(SGN*X0,SGN*YY) + (0.5 - 0.5*SGN)*PI
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
CCC      SIG0 = (SIG(JP) - SIG(JO))*DSIO
CCC      SIG1 = (SIG(JP) - SIG(JM))*DSIM
CCC      SSUM = SIG0 + SIG1
CCC      SDIF = SIG0 - SIG1
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
         QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF)
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
         QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF)
C
         DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*(DSIP+DSIO)
     &                                          - PDNI*(DSIP-DSIO))
         DQDM(JP) = DQDM(JP) + QOPI*( PSNI*DSIO - PDNI*DSIO)
         DQDM(JQ) = DQDM(JQ) + QOPI*( PSNI*DSIP + PDNI*DSIP)
C
        ENDIF
C
C------ calculate vortex panel contribution to Psi
        DXINV = 1.0/(X1-X2)
        PSIS = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
        PSID = ((X1+X2)*PSIS + 0.5*(RS2*G2-RS1*G1 + X1*X1-X2*X2))*DXINV
C
        PSX1 = 0.5*G1
        PSX2 = -.5*G2
        PSYY = T1-T2
C
        PDX1 = ((X1+X2)*PSX1 + PSIS - X1*G1 - PSID)*DXINV
        PDX2 = ((X1+X2)*PSX2 + PSIS + X2*G2 + PSID)*DXINV
        PDYY = ((X1+X2)*PSYY - YY*(G1-G2)         )*DXINV
C
        GSUM1 = GAMU(JP,1) + GAMU(JO,1)
        GSUM2 = GAMU(JP,2) + GAMU(JO,2)
        GDIF1 = GAMU(JP,1) - GAMU(JO,1)
        GDIF2 = GAMU(JP,2) - GAMU(JO,2)
C
        GSUM = GAM(JP) + GAM(JO)
        GDIF = GAM(JP) - GAM(JO)
C
        PSI = PSI + QOPI*(PSIS*GSUM + PSID*GDIF)
C
C------ dPsi/dGam
        DZDG(JO) = DZDG(JO) + QOPI*(PSIS-PSID)
        DZDG(JP) = DZDG(JP) + QOPI*(PSIS+PSID)
C
C------ dPsi/dni
        PSNI = PSX1*X1I + PSX2*X2I + PSYY*YYI
        PDNI = PDX1*X1I + PDX2*X2I + PDYY*YYI
        PSI_NI = PSI_NI + QOPI*(GSUM*PSNI + GDIF*PDNI)
C
        QTAN1 = QTAN1 + QOPI*(GSUM1*PSNI + GDIF1*PDNI)
        QTAN2 = QTAN2 + QOPI*(GSUM2*PSNI + GDIF2*PDNI)
C
        DQDG(JO) = DQDG(JO) + QOPI*(PSNI - PDNI)
        DQDG(JP) = DQDG(JP) + QOPI*(PSNI + PDNI)
C
        IF(GEOLIN) THEN
C
C------- dPsi/dn
         DZDN(JO) = DZDN(JO)+ QOPI*GSUM*(PSX1*X1O + PSX2*X2O + PSYY*YYO)
     &                      + QOPI*GDIF*(PDX1*X1O + PDX2*X2O + PDYY*YYO)
         DZDN(JP) = DZDN(JP)+ QOPI*GSUM*(PSX1*X1P + PSX2*X2P + PSYY*YYP)
     &                      + QOPI*GDIF*(PDX1*X1P + PDX2*X2P + PDYY*YYP)
C------- dPsi/dP
         Z_QDOF0 = Z_QDOF0
     &           + QOPI*((PSIS-PSID)*QF0(JO) + (PSIS+PSID)*QF0(JP))
         Z_QDOF1 = Z_QDOF1
     &           + QOPI*((PSIS-PSID)*QF1(JO) + (PSIS+PSID)*QF1(JP))
         Z_QDOF2 = Z_QDOF2
     &           + QOPI*((PSIS-PSID)*QF2(JO) + (PSIS+PSID)*QF2(JP))
         Z_QDOF3 = Z_QDOF3
     &           + QOPI*((PSIS-PSID)*QF3(JO) + (PSIS+PSID)*QF3(JP))
        ENDIF
C
C
   10 CONTINUE
C
   11 CONTINUE
      PSIG = 0.5*YY*(G1-G2) + X2*(T2-APAN) - X1*(T1-APAN)
      PGAM = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2)
C
      PSIGX1 = -(T1-APAN)
      PSIGX2 =   T2-APAN
      PSIGYY = 0.5*(G1-G2)
      PGAMX1 = 0.5*G1
      PGAMX2 = -.5*G2
      PGAMYY = T1-T2
C
      PSIGNI = PSIGX1*X1I + PSIGX2*X2I + PSIGYY*YYI
      PGAMNI = PGAMX1*X1I + PGAMX2*X2I + PGAMYY*YYI
C
C---- TE panel source and vortex strengths
      SIGTE1 = 0.5*SCS*(GAMU(JP,1) - GAMU(JO,1))
      SIGTE2 = 0.5*SCS*(GAMU(JP,2) - GAMU(JO,2))
      GAMTE1 = -.5*SDS*(GAMU(JP,1) - GAMU(JO,1))
      GAMTE2 = -.5*SDS*(GAMU(JP,2) - GAMU(JO,2))
C
      SIGTE = 0.5*SCS*(GAM(JP) - GAM(JO))
      GAMTE = -.5*SDS*(GAM(JP) - GAM(JO))
C
C---- TE panel contribution to Psi
      PSI = PSI + HOPI*(PSIG*SIGTE + PGAM*GAMTE)
C
C---- dPsi/dGam
      DZDG(JO) = DZDG(JO) - HOPI*PSIG*SCS*0.5
      DZDG(JP) = DZDG(JP) + HOPI*PSIG*SCS*0.5
C
      DZDG(JO) = DZDG(JO) + HOPI*PGAM*SDS*0.5
      DZDG(JP) = DZDG(JP) - HOPI*PGAM*SDS*0.5
C
C---- dPsi/dni
      PSI_NI = PSI_NI + HOPI*(PSIGNI*SIGTE + PGAMNI*GAMTE)
C
      QTAN1 = QTAN1 + HOPI*(PSIGNI*SIGTE1 + PGAMNI*GAMTE1)
      QTAN2 = QTAN2 + HOPI*(PSIGNI*SIGTE2 + PGAMNI*GAMTE2)
C
      DQDG(JO) = DQDG(JO) - HOPI*(PSIGNI*0.5*SCS - PGAMNI*0.5*SDS)
      DQDG(JP) = DQDG(JP) + HOPI*(PSIGNI*0.5*SCS - PGAMNI*0.5*SDS)
C
      IF(GEOLIN) THEN
C
C----- dPsi/dn
       DZDN(JO) = DZDN(JO)
     &          + HOPI*(PSIGX1*X1O + PSIGX2*X2O + PSIGYY*YYO)*SIGTE
     &          + HOPI*(PGAMX1*X1O + PGAMX2*X2O + PGAMYY*YYO)*GAMTE
       DZDN(JP) = DZDN(JP)
     &          + HOPI*(PSIGX1*X1P + PSIGX2*X2P + PSIGYY*YYP)*SIGTE
     &          + HOPI*(PGAMX1*X1P + PGAMX2*X2P + PGAMYY*YYP)*GAMTE
C
C----- dPsi/dP
       Z_QDOF0 = Z_QDOF0 + HOPI*PSIG*0.5*(QF0(JP)-QF0(JO))*SCS
     &                   - HOPI*PGAM*0.5*(QF0(JP)-QF0(JO))*SDS
       Z_QDOF1 = Z_QDOF1 + HOPI*PSIG*0.5*(QF1(JP)-QF1(JO))*SCS
     &                   - HOPI*PGAM*0.5*(QF1(JP)-QF1(JO))*SDS
       Z_QDOF2 = Z_QDOF2 + HOPI*PSIG*0.5*(QF2(JP)-QF2(JO))*SCS
     &                   - HOPI*PGAM*0.5*(QF2(JP)-QF2(JO))*SDS
       Z_QDOF3 = Z_QDOF3 + HOPI*PSIG*0.5*(QF3(JP)-QF3(JO))*SCS
     &                   - HOPI*PGAM*0.5*(QF3(JP)-QF3(JO))*SDS
C
      ENDIF
C
   12 CONTINUE
C
C**** Freestream terms
      PSI = PSI + QINF*(COSA*YI - SINA*XI)
C
C---- dPsi/dn
      PSI_NI = PSI_NI + QINF*(COSA*NYI - SINA*NXI)
C
      QTAN1 = QTAN1 + QINF*NYI
      QTAN2 = QTAN2 - QINF*NXI
C
C---- dPsi/dQinf
      Z_QINF = Z_QINF + (COSA*YI - SINA*XI)
C
C---- dPsi/dalfa
      Z_ALFA = Z_ALFA - QINF*(SINA*YI + COSA*XI)
C
C     DP note: In Xfoil there is more stuff below here IF(LIMAGE), 
C     but LIMAGE is always .FALSE., even in Xfoil
C
      RETURN
      END !PSILIN

C===================================================================70
C
C     Calculates two surface vorticity (gamma) distributions
C     for alpha = 0, 90  degrees.  These are superimposed
C     in SPECAL or SPECCL for specified alpha or CL.
C
C===================================================================70
      SUBROUTINE GGCALC

      use xfoil_inc
      use my_equivalence, only : my_equiv_3_2

C
C---- distance of internal control point ahead of sharp TE
C-    (fraction of smaller panel length adjacent to TE)
      BWT = 0.1
C
C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE)
     &  WRITE(*,*) 'Calculating unit vorticity distributions ...'
C
      DO 10 I=1, N
        GAM(I) = 0.
        GAMU(I,1) = 0.
        GAMU(I,2) = 0.
   10 CONTINUE
      PSIO = 0.
C
C---- Set up matrix system for  Psi = Psio  on airfoil surface.
C-    The unknowns are (dGamma)i and dPsio.
      DO 20 I=1, N
C
C------ calculate Psi and dPsi/dGamma array for current node
        CALL PSILIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_N,.FALSE.,.TRUE.)
C
        PSIINF = QINF*(COS(ALFA)*Y(I) - SIN(ALFA)*X(I))
C
C------ RES1 = PSI( 0) - PSIO
C------ RES2 = PSI(90) - PSIO
        RES1 =  QINF*Y(I)
        RES2 = -QINF*X(I)
C
C------ dRes/dGamma
        DO 201 J=1, N
          AIJ(I,J) = DZDG(J)
  201   CONTINUE
C
        DO 202 J=1, N
          BIJ(I,J) = -DZDM(J)

C         DP mod: copy to VM; used to replace equivalence statement
          call my_equiv_3_2(VM, BIJ, (/ 1, 1, 1 /), (/ 1, 1 /), 
     ?                     (/ I, J /), 1)

  202   CONTINUE
C
C------ dRes/dPsio
        AIJ(I,N+1) = -1.0
C
        GAMU(I,1) = -RES1
        GAMU(I,2) = -RES2
C
   20 CONTINUE
C
C---- set Kutta condition
C-    RES = GAM(1) + GAM(N)
      RES = 0.
C
      DO 30 J=1, N+1
        AIJ(N+1,J) = 0.0
   30 CONTINUE
C
      AIJ(N+1,1) = 1.0
      AIJ(N+1,N) = 1.0
C
      GAMU(N+1,1) = -RES
      GAMU(N+1,2) = -RES
C
C---- set up Kutta condition (no direct source influence)
      DO 32 J=1, N
        BIJ(N+1,J) = 0.

C       DP mod: copy to VM; used to replace equivalence statement
        call my_equiv_3_2(VM, BIJ, (/ 1, 1, 1 /), (/ 1, 1 /), 
     ?                   (/ N+1, J /), 1)

   32 CONTINUE
C
      IF(SHARP) THEN
C----- set zero internal velocity in TE corner 
C
C----- set TE bisector angle
       AG1 = ATAN2(-YP(1),-XP(1)    )
       AG2 = ATANC( YP(N), XP(N),AG1)
       ABIS = 0.5*(AG1+AG2)
       CBIS = COS(ABIS)
       SBIS = SIN(ABIS)
C
C----- minimum panel length adjacent to TE
       DS1 = SQRT( (X(1)-X(2)  )**2 + (Y(1)-Y(2)  )**2 )
       DS2 = SQRT( (X(N)-X(N-1))**2 + (Y(N)-Y(N-1))**2 )
       DSMIN = MIN( DS1 , DS2 )
C
C----- control point on bisector just ahead of TE point
       XBIS = XTE - BWT*DSMIN*CBIS
       YBIS = YTE - BWT*DSMIN*SBIS
ccc       write(*,*) xbis, ybis
C
C----- set velocity component along bisector line
       CALL PSILIN(0,XBIS,YBIS,-SBIS,CBIS,PSI,QBIS,.FALSE.,.TRUE.)
C
CCC--- RES = DQDGj*Gammaj + DQDMj*Massj + QINF*(COSA*CBIS + SINA*SBIS)
       RES = QBIS
C
C----- dRes/dGamma
       DO J=1, N
         AIJ(N,J) = DQDG(J)
       ENDDO
C
C----- -dRes/dMass
       DO J=1, N
         BIJ(N,J) = -DQDM(J)

C        DP mod: Copy to VM used to replace equivalence statement
         call my_equiv_3_2(VM, BIJ, (/ 1, 1, 1 /), (/ 1, 1 /), 
     ?                    (/ N, J /), 1)

       ENDDO
C
C----- dRes/dPsio
       AIJ(N,N+1) = 0.
C
C----- -dRes/dUinf
       GAMU(N,1) = -CBIS
C
C----- -dRes/dVinf
       GAMU(N,2) = -SBIS
C
      ENDIF
C
C---- LU-factor coefficient matrix AIJ
      CALL LUDCMP(IQX,N+1,AIJ,AIJPIV)
      LQAIJ = .TRUE.
C
C---- solve system for the two vorticity distributions
      CALL BAKSUB(IQX,N+1,AIJ,AIJPIV,GAMU(1,1))
      CALL BAKSUB(IQX,N+1,AIJ,AIJPIV,GAMU(1,2))
C
C---- set inviscid alpha=0,90 surface speeds for this geometry
      DO 50 I=1, N
        QINVU(I,1) = GAMU(I,1)
        QINVU(I,2) = GAMU(I,2)
   50 CONTINUE
C
      LGAMU = .TRUE.

      RETURN
      END

C===================================================================70
C
C     Sets inviscid panel tangential velocity for
C     current alpha.
C
C===================================================================70
      SUBROUTINE QISET

      use xfoil_inc
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
C
      DO 5 I=1, N+NW
        QINV  (I) =  COSA*QINVU(I,1) + SINA*QINVU(I,2)
        QINV_A(I) = -SINA*QINVU(I,1) + COSA*QINVU(I,2)
    5 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Integrates surface pressures to get CL and CM.
C     Integrates skin friction to get CDF.
C     Calculates dCL/dAlpha for prescribed-CL routines.
C
C===================================================================70
      SUBROUTINE CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, 
     &                  XREF,YREF,
     &                  CL,CM,CDP, CL_ALF,CL_MSQ)
      DIMENSION X(N),Y(N), GAM(N), GAM_A(N)
      REAL*8 MINF
C
ccC---- moment-reference coordinates
cc      XREF = 0.25
cc      YREF = 0.
C
      SA = SIN(ALFA)
      CA = COS(ALFA)
C
      BETA = SQRT(1.0 - MINF**2)
      BETA_MSQ = -0.5/BETA
C
      BFAC     = 0.5*MINF**2 / (1.0 + BETA)
      BFAC_MSQ = 0.5         / (1.0 + BETA)
     &         - BFAC        / (1.0 + BETA) * BETA_MSQ
C
      CL = 0.0
      CM = 0.0

      CDP = 0.0
C
      CL_ALF = 0.
      CL_MSQ = 0.
C
      I = 1
      CGINC = 1.0 - (GAM(I)/QINF)**2
      CPG1     = CGINC/(BETA + BFAC*CGINC)
      CPG1_MSQ = -CPG1/(BETA + BFAC*CGINC)*(BETA_MSQ + BFAC_MSQ*CGINC)
C
      CPI_GAM = -2.0*GAM(I)/QINF**2
      CPC_CPI = (1.0 - BFAC*CPG1)/ (BETA + BFAC*CGINC)
      CPG1_ALF = CPC_CPI*CPI_GAM*GAM_A(I)
C
      DO 10 I=1, N
        IP = I+1
        IF(I.EQ.N) IP = 1
C
        CGINC = 1.0 - (GAM(IP)/QINF)**2
        CPG2     = CGINC/(BETA + BFAC*CGINC)
        CPG2_MSQ = -CPG2/(BETA + BFAC*CGINC)*(BETA_MSQ + BFAC_MSQ*CGINC)
C
        CPI_GAM = -2.0*GAM(IP)/QINF**2
        CPC_CPI = (1.0 - BFAC*CPG2)/ (BETA + BFAC*CGINC)
        CPG2_ALF = CPC_CPI*CPI_GAM*GAM_A(IP)
C
        DX = (X(IP) - X(I))*CA + (Y(IP) - Y(I))*SA
        DY = (Y(IP) - Y(I))*CA - (X(IP) - X(I))*SA
        DG = CPG2 - CPG1
C
        AX = (0.5*(X(IP)+X(I))-XREF)*CA + (0.5*(Y(IP)+Y(I))-YREF)*SA
        AY = (0.5*(Y(IP)+Y(I))-YREF)*CA - (0.5*(X(IP)+X(I))-XREF)*SA
        AG = 0.5*(CPG2 + CPG1)
C
        DX_ALF = -(X(IP) - X(I))*SA + (Y(IP) - Y(I))*CA
        AG_ALF = 0.5*(CPG2_ALF + CPG1_ALF)
        AG_MSQ = 0.5*(CPG2_MSQ + CPG1_MSQ)
C
        CL     = CL     + DX* AG
        CDP    = CDP    - DY* AG
        CM     = CM     - DX*(AG*AX + DG*DX/12.0)
     &                  - DY*(AG*AY + DG*DY/12.0)
C
        CL_ALF = CL_ALF + DX*AG_ALF + AG*DX_ALF
        CL_MSQ = CL_MSQ + DX*AG_MSQ
C
        CPG1 = CPG2
        CPG1_ALF = CPG2_ALF
        CPG1_MSQ = CPG2_MSQ
   10 CONTINUE
C
      RETURN
      END ! CLCALC

C===================================================================70
C
C     Sets compressible Cp from speed.
C
C===================================================================70
      SUBROUTINE CPCALC(N,Q,QINF,MINF,CP,SILENT_MODE)

      DIMENSION Q(N),CP(N)
      REAL*8 MINF
      LOGICAL SILENT_MODE
C
      LOGICAL DENNEG
C
      BETA = SQRT(1.0 - MINF**2)
      BFAC = 0.5*MINF**2 / (1.0 + BETA)
C
      DENNEG = .FALSE.
C
      DO 20 I=1, N
        CPINC = 1.0 - (Q(I)/QINF)**2
        DEN = BETA + BFAC*CPINC
        CP(I) = CPINC / DEN
        IF(DEN .LE. 0.0) DENNEG = .TRUE.
  20  CONTINUE
C
C     DP mod: added SILENT_MODE option
      IF(DENNEG .AND. .NOT. SILENT_MODE) THEN
       WRITE(*,*)
       WRITE(*,*) 'CPCALC: Local speed too large. ',
     &            'Compressibility corrections invalid.'
      ENDIF
C
      RETURN
      END ! CPCALC

C===================================================================70
C
C Sets MINF from input
C
C===================================================================70
      SUBROUTINE MINFSET(MACH_INPUT)

      use xfoil_inc

      REAL*8 MACH_INPUT

      IF(MINF1.GE.1.0) THEN
        WRITE(*,*) 'Supersonic freestream not allowed'
        STOP
      ENDIF
      MINF1 = MACH_INPUT

      CALL MRCL(1.0,MINF_CL,REINF_CL)
      CALL COMSET

      IF ( (MINF.GT.0.0) .AND. (.NOT. SILENT_MODE) ) 
     &   WRITE(*,1300) CPSTAR, QSTAR/QINF
 1300 FORMAT(/' Sonic Cp =', F10.2, '      Sonic Q/Qinf =', F10.3/)

      CALL CPCALC(N,QINV,QINF,MINF,CPI,SILENT_MODE)
      IF(LVISC) CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV,SILENT_MODE)
      LVCONV = .FALSE.

      END ! MINFSET

C===================================================================70
C
C     Converges to specified alpha.
C
C===================================================================70
      SUBROUTINE SPECAL

      use xfoil_inc

      REAL*8 MINF_CLM, MSQ_CLM

C
C---- calculate surface vorticity distributions for alpha = 0, 90 degrees
      IF(.NOT.LGAMU .OR. .NOT.LQAIJ) CALL GGCALC
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
C
C---- superimpose suitably weighted  alpha = 0, 90  distributions
      DO 50 I=1, N
        GAM(I)   =  COSA*GAMU(I,1) + SINA*GAMU(I,2)
        GAM_A(I) = -SINA*GAMU(I,1) + COSA*GAMU(I,2)
   50 CONTINUE
      PSIO = COSA*GAMU(N+1,1) + SINA*GAMU(N+1,2)
C
      CALL TECALC
      CALL QISET
C
C---- set initial guess for the Newton variable CLM
      CLM = 1.0
C
C---- set corresponding  M(CLM), Re(CLM)
      CALL MRCL(CLM,MINF_CLM,REINF_CLM)
      CALL COMSET
C
C---- set corresponding CL(M)
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &            CL,CM,CDP, CL_ALF,CL_MSQ)
C
C---- iterate on CLM
      DO 100 ITCL=1, 20
C
        MSQ_CLM = 2.0*MINF*MINF_CLM
        DCLM = (CL - CLM)/(1.0 - CL_MSQ*MSQ_CLM)
C
        CLM1 = CLM
        RLX = 1.0
C
C------ under-relaxation loop to avoid driving M(CL) above 1
        DO 90 IRLX=1, 12
C
          CLM = CLM1 + RLX*DCLM
C
C-------- set new freestream Mach M(CLM)
          CALL MRCL(CLM,MINF_CLM,REINF_CLM)
C
C-------- if Mach is OK, go do next Newton iteration
          IF(MATYP.EQ.1 .OR. MINF.EQ.0.0 .OR. MINF_CLM.NE.0.0) GO TO 91
C
          RLX = 0.5*RLX
   90   CONTINUE
   91   CONTINUE
C
C------ set new CL(M)
        CALL COMSET
        CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &              CL,CM,CDP,CL_ALF,CL_MSQ)
C
        IF(ABS(DCLM).LE.1.0E-6) GO TO 110
C
  100 CONTINUE
C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE) 
     &  WRITE(*,*) 'SPECAL:  Minf convergence failed'
  110 CONTINUE
C
C---- set final Mach, CL, Cp distributions, and hinge moment
      CALL MRCL(CL,MINF_CL,REINF_CL)
      CALL COMSET
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &            CL,CM,CDP, CL_ALF,CL_MSQ)
      CALL CPCALC(N,QINV,QINF,MINF,CPI,SILENT_MODE)
      IF(LVISC) THEN
       CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV,SILENT_MODE)
       CALL CPCALC(N+NW,QINV,QINF,MINF,CPI,SILENT_MODE)
      ELSE
       CALL CPCALC(N,QINV,QINF,MINF,CPI,SILENT_MODE)
      ENDIF
C      IF(LFLAP) CALL MHINGE
C
      RETURN
      END ! SPECAL

C===================================================================70
C
C     Converges to specified inviscid CL.
C
C===================================================================70
      SUBROUTINE SPECCL

      use xfoil_inc
C
C---- calculate surface vorticity distributions for alpha = 0, 90 degrees
      IF(.NOT.LGAMU .OR. .NOT.LQAIJ) CALL GGCALC
C
C---- set freestream Mach from specified CL -- Mach will be held fixed
      CALL MRCL(CLSPEC,MINF_CL,REINF_CL)
      CALL COMSET
C
C---- current alpha is the initial guess for Newton variable ALFA
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
      DO 10 I=1, N
        GAM(I)   =  COSA*GAMU(I,1) + SINA*GAMU(I,2)
        GAM_A(I) = -SINA*GAMU(I,1) + COSA*GAMU(I,2)
   10 CONTINUE
      PSIO = COSA*GAMU(N+1,1) + SINA*GAMU(N+1,2)
C
C---- get corresponding CL, CL_alpha, CL_Mach
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &            CL,CM,CDP, CL_ALF,CL_MSQ)
C
C---- Newton loop for alpha to get specified inviscid CL
      DO 100 ITAL=1, 20
C
        DALFA = (CLSPEC - CL) / CL_ALF
        RLX = 1.0
C
        ALFA = ALFA + RLX*DALFA
C
C------ set new surface speed distribution
        COSA = COS(ALFA)
        SINA = SIN(ALFA)
        DO 40 I=1, N
          GAM(I)   =  COSA*GAMU(I,1) + SINA*GAMU(I,2)
          GAM_A(I) = -SINA*GAMU(I,1) + COSA*GAMU(I,2)
   40   CONTINUE
        PSIO = COSA*GAMU(N+1,1) + SINA*GAMU(N+1,2)
C
C------ set new CL(alpha)
        CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &              CL,CM,CDP,CL_ALF,CL_MSQ)
C
        IF(ABS(DALFA).LE.1.0E-6) GO TO 110
  100 CONTINUE
C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE) 
     &  WRITE(*,*) 'SPECCL:  CL convergence failed'
  110 CONTINUE
C
C---- set final surface speed and Cp distributions
      CALL TECALC
      CALL QISET
      IF(LVISC) THEN
       CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV,SILENT_MODE)
       CALL CPCALC(N+NW,QINV,QINF,MINF,CPI,SILENT_MODE)
      ELSE
       CALL CPCALC(N,QINV,QINF,MINF,CPI,SILENT_MODE)
      ENDIF
C      IF(LFLAP) CALL MHINGE
C
      RETURN
      END ! SPECCL

C===================================================================70
C
C     Sets inviscid tangential velocity for alpha = 0, 90
C     on wake due to freestream and airfoil surface vorticity.
C
C===================================================================70
      SUBROUTINE QWCALC

      use xfoil_inc
C
C---- first wake point (same as TE)
      QINVU(N+1,1) = QINVU(N,1)
      QINVU(N+1,2) = QINVU(N,2)
C
C---- rest of wake
      DO 10 I=N+2, N+NW
        CALL PSILIN(I,X(I),Y(I),NX(I),NY(I),PSI,PSI_NI,.FALSE.,.FALSE.)
        QINVU(I,1) = QTAN1
        QINVU(I,2) = QTAN2
   10 CONTINUE
C
      RETURN
      END

C===================================================================70
C===================================================================70
      SUBROUTINE GAMQV
      
      use xfoil_inc
C
      DO 10 I=1, N
        GAM(I)   = QVIS(I)
        GAM_A(I) = QINV_A(I)
   10 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Locates stagnation point arc length 
C     location SST and panel index IST.
C
C===================================================================70
      SUBROUTINE STFIND

      use xfoil_inc
C
      DO 10 I=1, N-1
        IF(GAM(I).GE.0.0 .AND. GAM(I+1).LT.0.0) GO TO 11
   10 CONTINUE
C
C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE)
     &  WRITE(*,*) 'STFIND: Stagnation point not found. Continuing ...'
      I = N/2
C
   11 CONTINUE
C
      IST = I
      DGAM = GAM(I+1) - GAM(I)
      DS = S(I+1) - S(I)
C
C---- evaluate so as to minimize roundoff for very small GAM(I) or GAM(I+1)
      IF(GAM(I) .LT. -GAM(I+1)) THEN
       SST = S(I)   - DS*(GAM(I)  /DGAM)
      ELSE
       SST = S(I+1) - DS*(GAM(I+1)/DGAM)
      ENDIF
C
C---- tweak stagnation point if it falls right on a node (very unlikely)
      IF(SST .LE. S(I)  ) SST = S(I)   + 1.0E-7
      IF(SST .GE. S(I+1)) SST = S(I+1) - 1.0E-7
C
      SST_GO = (SST  - S(I+1))/DGAM
      SST_GP = (S(I) - SST   )/DGAM
C
      RETURN
      END

C===================================================================70
C
C     Sets  BL location -> panel location  pointer array IPAN
C
C===================================================================70
      SUBROUTINE IBLPAN

      use xfoil_inc
C
C---- top surface first
      IS = 1
C
      IBL = 1
      DO 10 I=IST, 1, -1
        IBL = IBL+1
        IPAN(IBL,IS) = I
        VTI(IBL,IS) = 1.0
   10 CONTINUE
C
      IBLTE(IS) = IBL
      NBL(IS) = IBL
C
C---- bottom surface next
      IS = 2
C
      IBL = 1
      DO 20 I=IST+1, N
        IBL = IBL+1
        IPAN(IBL,IS) = I
        VTI(IBL,IS) = -1.0
   20 CONTINUE
C
C---- wake
      IBLTE(IS) = IBL
C
      DO 25 IW=1, NW
        I = N+IW
        IBL = IBLTE(IS)+IW
        IPAN(IBL,IS) = I
         VTI(IBL,IS) = -1.0
   25 CONTINUE
C
      NBL(IS) = IBLTE(IS) + NW
C     DP note: stopped here
C
C---- upper wake pointers (for plotting only)
      DO 35 IW=1, NW
        IPAN(IBLTE(1)+IW,1) = IPAN(IBLTE(2)+IW,2)
         VTI(IBLTE(1)+IW,1) = 1.0
   35 CONTINUE
C
C
      IBLMAX = MAX(IBLTE(1),IBLTE(2)) + NW
      IF(IBLMAX.GT.IVX) THEN
        WRITE(*,*) ' ***  BL array overflow.'
        WRITE(*,*) ' ***  Increase IVX to at least', IBLMAX
        STOP
      ENDIF
C
      LIPAN = .TRUE.
      RETURN
      END

C===================================================================70
C
C     Sets BL arc length array on each airfoil side and wake
C
C===================================================================70
      SUBROUTINE XICALC

      use xfoil_inc

      DATA XFEPS / 1.0E-7 /
C
C---- minimum xi node arc length near stagnation point
      XEPS = XFEPS*(S(N)-S(1))
C
      IS = 1
C
      XSSI(1,IS) = 0.
C
      DO 10 IBL=2, IBLTE(IS)
        I = IPAN(IBL,IS)
        XSSI(IBL,IS) = MAX( SST - S(I) , XEPS )
   10 CONTINUE
C
C
      IS = 2
C
      XSSI(1,IS) = 0.
C
      DO 20 IBL=2, IBLTE(IS)
        I = IPAN(IBL,IS)
        XSSI(IBL,IS) = MAX( S(I) - SST , XEPS )
   20 CONTINUE
C
C
      IS1 = 1
      IS2 = 2
C
      IBL1 = IBLTE(IS1) + 1
      XSSI(IBL1,IS1) = XSSI(IBL1-1,IS1)
C
      IBL2 = IBLTE(IS2) + 1
      XSSI(IBL2,IS2) = XSSI(IBL2-1,IS2)
C
      DO 25 IBL=IBLTE(IS)+2, NBL(IS)
        I = IPAN(IBL,IS)
        DXSSI = SQRT((X(I)-X(I-1))**2 + (Y(I)-Y(I-1))**2)
C
        IBL1 = IBLTE(IS1) + IBL - IBLTE(IS)
        IBL2 = IBLTE(IS2) + IBL - IBLTE(IS)
        XSSI(IBL1,IS1) = XSSI(IBL1-1,IS1) + DXSSI
        XSSI(IBL2,IS2) = XSSI(IBL2-1,IS2) + DXSSI
   25 CONTINUE
C
C---- trailing edge flap length to TE gap ratio
      TELRAT = 2.50
C
C---- set up parameters for TE flap cubics
C
ccc   DWDXTE = YP(1)/XP(1) + YP(N)/XP(N)    !!! BUG  2/2/95
C
      CROSP = (XP(1)*YP(N) - YP(1)*XP(N))
     &      / SQRT(  (XP(1)**2 + YP(1)**2)
     &              *(XP(N)**2 + YP(N)**2) )
      DWDXTE = CROSP / SQRT(1.0 - CROSP**2)
C
C---- limit cubic to avoid absurd TE gap widths
      DWDXTE = MAX(DWDXTE,-3.0/TELRAT)
      DWDXTE = MIN(DWDXTE, 3.0/TELRAT)
C
      AA =  3.0 + TELRAT*DWDXTE
      BB = -2.0 - TELRAT*DWDXTE
C
      IF(SHARP) THEN
       DO 30 IW=1, NW
         WGAP(IW) = 0.
   30  CONTINUE
      ELSE
C----- set TE flap (wake gap) array
       IS = 2
       DO 35 IW=1, NW
         IBL = IBLTE(IS) + IW
         ZN = 1.0 - (XSSI(IBL,IS)-XSSI(IBLTE(IS),IS)) / (TELRAT*ANTE)
         WGAP(IW) = 0.
         IF(ZN.GE.0.0) WGAP(IW) = ANTE * (AA + BB*ZN)*ZN**2
   35  CONTINUE
      ENDIF
C
      RETURN
      END

C===================================================================70
C
C     Sets inviscid Ue from panel inviscid tangential velocity
C
C===================================================================70
      SUBROUTINE UICALC

      use xfoil_inc
C
      DO 10 IS=1, 2
        UINV  (1,IS) = 0.
        UINV_A(1,IS) = 0.
        DO 110 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
          UINV  (IBL,IS) = VTI(IBL,IS)*QINV  (I)
          UINV_A(IBL,IS) = VTI(IBL,IS)*QINV_A(I)
  110   CONTINUE
   10 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Sets panel viscous tangential velocity from viscous Ue
C
C===================================================================70
      SUBROUTINE QVFUE

      use xfoil_inc
C
      DO 1 IS=1, 2
        DO 10 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
          QVIS(I) = VTI(IBL,IS)*UEDG(IBL,IS)
   10   CONTINUE
    1 CONTINUE
C
      RETURN
      END

C===================================================================70
C===================================================================70
      SUBROUTINE CDCALC

      use xfoil_inc
C
      SA = SIN(ALFA)
      CA = COS(ALFA)
C
      IF(LVISC .AND. LBLINI) THEN
C
C----- set variables at the end of the wake
       THWAKE = THET(NBL(2),2)
       URAT   = UEDG(NBL(2),2)/QINF
       UEWAKE = UEDG(NBL(2),2) * (1.0-TKLAM) / (1.0 - TKLAM*URAT**2)
       SHWAKE = DSTR(NBL(2),2)/THET(NBL(2),2)
C
C----- extrapolate wake to downstream infinity using Squire-Young relation
C      (reduces errors of the wake not being long enough)
       CD = 2.0*THWAKE * (UEWAKE/QINF)**(0.5*(5.0+SHWAKE))
C
      ELSE
C
       CD = 0.0
C
      ENDIF
C
C---- calculate friction drag coefficient
      CDF = 0.0
      DO 20 IS=1, 2
        DO 205 IBL=3, IBLTE(IS)
          I  = IPAN(IBL  ,IS)
          IM = IPAN(IBL-1,IS)
          DX = (X(I) - X(IM))*CA + (Y(I) - Y(IM))*SA
          CDF = CDF + 0.5*(TAU(IBL,IS)+TAU(IBL-1,IS))*DX * 2.0/QINF**2
 205    CONTINUE
 20   CONTINUE
C
      RETURN
      END ! CDCALC

C===================================================================70
C
C     Moves stagnation point location to new panel.
C
C===================================================================70
      SUBROUTINE STMOVE

      use xfoil_inc
C
C---- locate new stagnation point arc length SST from GAM distribution
      ISTOLD = IST
      CALL STFIND
C
      IF(ISTOLD.EQ.IST) THEN
C
C----- recalculate new arc length array
       CALL XICALC
C
      ELSE
C
CCC       WRITE(*,*) 'STMOVE: Resetting stagnation point'
C
C----- set new BL position -> panel position  pointers
       CALL IBLPAN
C
C----- set new inviscid BL edge velocity UINV from QINV
       CALL UICALC
C
C----- recalculate new arc length array
       CALL XICALC
C
C----- set  BL position -> system line  pointers
       CALL IBLSYS
C
       IF(IST.GT.ISTOLD) THEN
C------ increase in number of points on top side (IS=1)
        IDIF = IST-ISTOLD
C
        ITRAN(1) = ITRAN(1) + IDIF
        ITRAN(2) = ITRAN(2) - IDIF
C
C------ move top side BL variables downstream
        DO 110 IBL=NBL(1), IDIF+2, -1
          CTAU(IBL,1) = CTAU(IBL-IDIF,1)
          THET(IBL,1) = THET(IBL-IDIF,1)
          DSTR(IBL,1) = DSTR(IBL-IDIF,1)
          UEDG(IBL,1) = UEDG(IBL-IDIF,1)
  110   CONTINUE            
C
C------ set BL variables between old and new stagnation point
        DUDX = UEDG(IDIF+2,1)/XSSI(IDIF+2,1)
        DO 115 IBL=IDIF+1, 2, -1
          CTAU(IBL,1) = CTAU(IDIF+2,1)
          THET(IBL,1) = THET(IDIF+2,1)
          DSTR(IBL,1) = DSTR(IDIF+2,1)
          UEDG(IBL,1) = DUDX * XSSI(IBL,1)
  115   CONTINUE
C
C------ move bottom side BL variables upstream
        DO 120 IBL=2, NBL(2)
          CTAU(IBL,2) = CTAU(IBL+IDIF,2)
          THET(IBL,2) = THET(IBL+IDIF,2)
          DSTR(IBL,2) = DSTR(IBL+IDIF,2)
          UEDG(IBL,2) = UEDG(IBL+IDIF,2)
  120   CONTINUE            
C
       ELSE
C------ increase in number of points on bottom side (IS=2)
        IDIF = ISTOLD-IST
C
        ITRAN(1) = ITRAN(1) - IDIF
        ITRAN(2) = ITRAN(2) + IDIF
C
C------ move bottom side BL variables downstream
        DO 210 IBL=NBL(2), IDIF+2, -1
          CTAU(IBL,2) = CTAU(IBL-IDIF,2)
          THET(IBL,2) = THET(IBL-IDIF,2)
          DSTR(IBL,2) = DSTR(IBL-IDIF,2)
          UEDG(IBL,2) = UEDG(IBL-IDIF,2)
  210   CONTINUE            
C
C------ set BL variables between old and new stagnation point
        DUDX = UEDG(IDIF+2,2)/XSSI(IDIF+2,2)


c        write(*,*) 'idif Ue xi dudx', 
c     &    idif, UEDG(idif+2,2), xssi(idif+2,2), dudx

        DO 215 IBL=IDIF+1, 2, -1
          CTAU(IBL,2) = CTAU(IDIF+2,2)
          THET(IBL,2) = THET(IDIF+2,2)
          DSTR(IBL,2) = DSTR(IDIF+2,2)
          UEDG(IBL,2) = DUDX * XSSI(IBL,2)
  215   CONTINUE

c        write(*,*) 'Uenew xinew', idif+1, uedg(idif+1,2), xssi(idif+1,2)

C
C------ move top side BL variables upstream
        DO 220 IBL=2, NBL(1)
          CTAU(IBL,1) = CTAU(IBL+IDIF,1)
          THET(IBL,1) = THET(IBL+IDIF,1)
          DSTR(IBL,1) = DSTR(IBL+IDIF,1)
          UEDG(IBL,1) = UEDG(IBL+IDIF,1)
  220   CONTINUE            
       ENDIF
C
C----- tweak Ue so it's not zero, in case stag. point is right on node
       UEPS = 1.0E-7
       DO IS = 1, 2
         DO IBL = 2, NBL(IS)
           I = IPAN(IBL,IS)
           IF(UEDG(IBL,IS).LE.UEPS) THEN
            UEDG(IBL,IS) = UEPS
            QVIS(I) = VTI(IBL,IS)*UEPS
            GAM(I)  = VTI(IBL,IS)*UEPS
           ENDIF
         ENDDO
       ENDDO
C
      ENDIF
C
C---- set new mass array since Ue has been tweaked
      DO 50 IS=1, 2
        DO 510 IBL=2, NBL(IS)
          MASS(IBL,IS) = DSTR(IBL,IS)*UEDG(IBL,IS)
  510   CONTINUE
   50 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Converges viscous operating point
C
C===================================================================70
      SUBROUTINE VISCAL(NITER1)

      use xfoil_inc
C
C---- convergence tolerance
      DATA EPS1 / 1.0E-4 /
C
      NITER = NITER1

C     DP mod: variable to notify of infinite loop condition (and halt)
      XFOIL_FAIL = .FALSE.
C
C---- calculate wake trajectory from current inviscid solution if necessary
      IF(.NOT.LWAKE) THEN
       CALL XYWAKE
      ENDIF
C
C---- set velocities on wake from airfoil vorticity for alpha=0, 90
      CALL QWCALC
C
C---- set velocities on airfoil and wake for initial alpha
      CALL QISET
C
      IF(.NOT.LIPAN) THEN
C
       IF(LBLINI) CALL GAMQV
C
C----- locate stagnation point arc length position and panel index
       CALL STFIND
C
C----- set  BL position -> panel position  pointers
       CALL IBLPAN
C
C----- calculate surface arc length array for current stagnation point location
       CALL XICALC
C
C----- set  BL position -> system line  pointers
       CALL IBLSYS
C
      ENDIF
C
C---- set inviscid BL edge velocity UINV from QINV
      CALL UICALC
C
      IF(.NOT.LBLINI) THEN
C
C----- set initial Ue from inviscid Ue
       DO IBL=1, NBL(1)
         UEDG(IBL,1) = UINV(IBL,1)
       ENDDO
C
       DO IBL=1, NBL(2)
         UEDG(IBL,2) = UINV(IBL,2)
       ENDDO
C
      ENDIF
C
      IF(LVCONV) THEN
C----- set correct CL if converged point exists
       CALL QVFUE
       IF(LVISC) THEN
        CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV,SILENT_MODE)
        CALL CPCALC(N+NW,QINV,QINF,MINF,CPI,SILENT_MODE)
       ELSE
        CALL CPCALC(N,QINV,QINF,MINF,CPI,SILENT_MODE)
       ENDIF
       CALL GAMQV
       CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &             CL,CM,CDP, CL_ALF,CL_MSQ)
       CALL CDCALC
      ENDIF
C
C---- set up source influence matrix if it doesn't exist
      IF(.NOT.LWDIJ .OR. .NOT.LADIJ) CALL QDCALC
C
C---- Newton iteration for entire BL solution
C     DP mod: set default NITER to 10
C      IF(NITER.EQ.0) CALL ASKI('Enter number of iterations^',NITER)
      IF(NITER.EQ.0) NITER = 10
C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE) THEN
        WRITE(*,*)
        WRITE(*,*) 'Solving BL system ...'
      END IF
      DO 1000 ITER=1, NITER
C
C------ fill Newton system for BL variables
        CALL SETBL
C       DP mod: check for infinite loop condition
        IF (XFOIL_FAIL) THEN
          CL = -0.1
          CD = 1000.0
          CM = -10.0
          RMSBL = 1000.0
          RETURN
        ENDIF
C
C------ solve Newton system with custom solver
        CALL BLSOLV
C
C------ update BL variables
        CALL UPDATE
C
        IF(LALFA) THEN
C------- set new freestream Mach, Re from new CL
         CALL MRCL(CL,MINF_CL,REINF_CL)
         CALL COMSET
        ELSE
C------- set new inviscid speeds QINV and UINV for new alpha
         CALL QISET
         CALL UICALC
        ENDIF
C
C------ calculate edge velocities QVIS(.) from UEDG(..)
        CALL QVFUE
C
C------ set GAM distribution from QVIS
        CALL GAMQV
C
C------ relocate stagnation point
        CALL STMOVE
C
C------ set updated CL,CD
        CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XCMREF,YCMREF,
     &              CL,CM,CDP,CL_ALF,CL_MSQ)
        CALL CDCALC
C
C------ display changes and test for convergence
C       DP mod: added SILENT_MODE options
        IF(.NOT. SILENT_MODE .AND. RLX.LT.1.0) 
     &   WRITE(*,2000) ITER, RMSBL, RMXBL, VMXBL,IMXBL,ISMXBL,RLX
        IF(.NOT. SILENT_MODE .AND. RLX.EQ.1.0) 
     &   WRITE(*,2010) ITER, RMSBL, RMXBL, VMXBL,IMXBL,ISMXBL
         CDPDIF = CD - CDF
         IF(.NOT. SILENT_MODE)
     &     WRITE(*,2020) ALFA/DTOR, CL, CM, CD, CDF, CDPDIF
C
        IF(RMSBL .LT. EPS1) THEN
         LVCONV = .TRUE.
         AVISC = ALFA
         MVISC = MINF
         GO TO 90
        ENDIF
C
 1000 CONTINUE
C     DP mod: added SILENT_MODE option
      IF(.NOT. SILENT_MODE) WRITE(*,*) 'VISCAL:  Convergence failed'
C
   90 CONTINUE
      CALL CPCALC(N+NW,QINV,QINF,MINF,CPI,SILENT_MODE)
      CALL CPCALC(N+NW,QVIS,QINF,MINF,CPV,SILENT_MODE)
C      IF(LFLAP) CALL MHINGE
      RETURN
C....................................................................
 2000   FORMAT
     &   (/1X,I3,'   rms: ',E10.4,'   max: ',E10.4,3X,A1,' at ',I4,I3,
     &     '   RLX:',F6.3)
 2010   FORMAT
     &   (/1X,I3,'   rms: ',E10.4,'   max: ',E10.4,3X,A1,' at ',I4,I3)
 2020   FORMAT
     &   ( 1X,3X,'   a =', F7.3,'      CL =',F8.4  /
     &     1X,3X,'  Cm =', F8.4, '     CD =',F9.5,
     &           '   =>   CDf =',F9.5,'    CDp =',F9.5)
      END ! VISCAL
