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
C
C     Sets actual Mach, Reynolds numbers
C     from unit-CL values and specified CLS
C     depending on MATYP,RETYP flags.
C
C===================================================================70
      SUBROUTINE MRCL(CLS,M_CLS,R_CLS)

      use xfoil_inc

      REAL*8 M_CLS
C
      CLA = MAX( CLS , 0.000001 )
C
C     DP mod: added SILENT_MODE option
      IF(RETYP.LT.1 .OR. RETYP.GT.3) THEN
        IF (.NOT. SILENT_MODE) THEN
          WRITE(*,*) 'MRCL:  Illegal Re(CL) dependence trigger.'
          WRITE(*,*) '       Setting fixed Re.'
        ENDIF
        RETYP = 1
      ENDIF
      IF(MATYP.LT.1 .OR. MATYP.GT.3) THEN
        IF (.NOT. SILENT_MODE) THEN
          WRITE(*,*) 'MRCL:  Illegal Mach(CL) dependence trigger.'
          WRITE(*,*) '       Setting fixed Mach.'
        ENDIF
        MATYP = 1
      ENDIF
C
C
      IF(MATYP.EQ.1) THEN
C
        MINF  = MINF1
        M_CLS = 0.
C
      ELSE IF(MATYP.EQ.2) THEN
C
        MINF  =  MINF1/SQRT(CLA)
        M_CLS = -0.5*MINF/CLA
C
      ELSE IF(MATYP.EQ.3) THEN
C
        MINF  = MINF1
        M_CLS = 0.
C
      ENDIF
C
C
      IF(RETYP.EQ.1) THEN
C
        REINF = REINF1
        R_CLS = 0.
C
      ELSE IF(RETYP.EQ.2) THEN
C
        REINF =  REINF1/SQRT(CLA)
        R_CLS = -0.5*REINF/CLA
C
      ELSE IF(RETYP.EQ.3) THEN
C
        REINF =  REINF1/CLA
        R_CLS = -REINF /CLA
C
      ENDIF
C
C
      IF(MINF .GE. 0.99) THEN
        IF (.NOT. SILENT_MODE) THEN
          WRITE(*,*)
          WRITE(*,*) 'MRCL: CL too low for chosen Mach(CL) dependence'
          WRITE(*,*) '      Aritificially limiting Mach to  0.99'
        ENDIF
        MINF = 0.99
        M_CLS = 0.
      ENDIF
C
      RRAT = 1.0
      IF(REINF1 .GT. 0.0) RRAT = REINF/REINF1
C
      IF(RRAT .GT. 100.0) THEN
        IF (.NOT. SILENT_MODE) THEN
          WRITE(*,*)
          WRITE(*,*) 'MRCL: CL too low for chosen Re(CL) dependence'
          WRITE(*,*) '      Aritificially limiting Re to ',REINF1*100.0
        ENDIF
        REINF = REINF1*100.0
        R_CLS = 0.
      ENDIF
C
      RETURN
      END ! MRCL

C===================================================================70
C
C     Sets the BL Newton system line number
C     corresponding to each BL station.
C
C===================================================================70
      SUBROUTINE IBLSYS

      use xfoil_inc
      use blpar_inc
      use xbl_inc
C
      IV = 0
      DO 10 IS=1, 2
        DO 110 IBL=2, NBL(IS)
          IV = IV+1
          ISYS(IBL,IS) = IV
  110   CONTINUE
   10 CONTINUE
C
      NSYS = IV
      IF(NSYS.GT.2*IVX) STOP '*** IBLSYS: BL system array overflow. ***'
C
      RETURN
      END

C===================================================================70
C 
C Set Karman-Tsien parameter TKLAM
C
C===================================================================70
      SUBROUTINE COMSET

      use xfoil_inc
C
      BETA = SQRT(1.0 - MINF**2)
      BETA_MSQ = -0.5/BETA
C
      TKLAM   = MINF**2 / (1.0 + BETA)**2
      TKL_MSQ =     1.0 / (1.0 + BETA)**2
     &    - 2.0*TKLAM/ (1.0 + BETA) * BETA_MSQ
C
C---- set sonic Pressure coefficient and speed
      IF(MINF.EQ.0.0) THEN
       CPSTAR = -999.0
       QSTAR = 999.0
      ELSE
       CPSTAR = 2.0 / (GAMMA*MINF**2)
     &        * (( (1.0 + 0.5*GAMM1*MINF**2)
     &            /(1.0 + 0.5*GAMM1        ))**(GAMMA/GAMM1) - 1.0)
       QSTAR = QINF/MINF
     &       * SQRT( (1.0 + 0.5*GAMM1*MINF**2)
     &              /(1.0 + 0.5*GAMM1        ) )
      ENDIF
C
      RETURN
      END ! COMSET

C===================================================================70
C
C     Sets forced-transition BL coordinate locations.
C
C===================================================================70
      SUBROUTINE XIFSET(IS)

      use xfoil_inc
      use blpar_inc
      use xbl_inc
C
      IF(XSTRIP(IS).GE.1.0) THEN
       XIFORC = XSSI(IBLTE(IS),IS)
       RETURN
      ENDIF
C
      CHX = XTE - XLE
      CHY = YTE - YLE
      CHSQ = CHX**2 + CHY**2
C
C---- calculate chord-based x/c, y/c
      DO 10 I=1, N
        W1(I) = ((X(I)-XLE)*CHX + (Y(I)-YLE)*CHY) / CHSQ
        W2(I) = ((Y(I)-YLE)*CHX - (X(I)-XLE)*CHY) / CHSQ
 10   CONTINUE
C
      CALL SPLIND(W1,W3,S,N,-999.0,-999.0)
      CALL SPLIND(W2,W4,S,N,-999.0,-999.0)
C
      IF(IS.EQ.1) THEN
C
C----- set approximate arc length of forced transition point for SINVRT
       STR = SLE + (S(1)-SLE)*XSTRIP(IS)
C
C----- calculate actual arc length
       CALL SINVRT(STR,XSTRIP(IS),W1,W3,S,N,SILENT_MODE)
C
C----- set BL coordinate value
       XIFORC = MIN( (SST - STR) , XSSI(IBLTE(IS),IS) )
C
      ELSE
C----- same for bottom side
C
       STR = SLE + (S(N)-SLE)*XSTRIP(IS)
       CALL SINVRT(STR,XSTRIP(IS),W1,W3,S,N,SILENT_MODE)
       XIFORC = MIN( (STR - SST) , XSSI(IBLTE(IS),IS) )
C
      ENDIF
C
      IF(XIFORC .LT. 0.0) THEN
C      DP mod: added SILENT_MODE option
       IF (.NOT. SILENT_MODE) WRITE(*,1000) IS
 1000  FORMAT(/' ***  Stagnation point is past trip on side',I2,'  ***')
       XIFORC = XSSI(IBLTE(IS),IS)
      ENDIF
C
      RETURN
      END

C===================================================================70
C
C     Set BL primary "2" variables from parameter list
C
C===================================================================70
      SUBROUTINE BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)

      use xbl_inc
      IMPLICIT REAL(M)
C
      X2 = XSI
      AMPL2 = AMI
      S2  = CTI
      T2  = THI
      D2  = DSI - DSWAKI
      DW2 = DSWAKI
C
      U2 = UEI*(1.0-TKBL) / (1.0 - TKBL*(UEI/QINFBL)**2)
      U2_UEI = (1.0 + TKBL*(2.0*U2*UEI/QINFBL**2 - 1.0))
     &       / (1.0 - TKBL*(UEI/QINFBL)**2)
      U2_MS  = (U2*(UEI/QINFBL)**2  -  UEI)*TKBL_MS
     &                    / (1.0 - TKBL*(UEI/QINFBL)**2)
C
      RETURN
      END ! BLPRV

C===================================================================70
C===================================================================70
      SUBROUTINE HKIN( H, MSQ, HK, HK_H, HK_MSQ )

      REAL*8 :: MSQ
C
C---- calculate kinematic shape parameter (assuming air)
C     (from Whitfield )
      HK     =    (H - 0.29*MSQ)/(1.0 + 0.113*MSQ)
      HK_H   =     1.0          /(1.0 + 0.113*MSQ)
      HK_MSQ = (-.29 - 0.113*HK)/(1.0 + 0.113*MSQ)
C
      RETURN
      END
C===================================================================70
C
C     Calculates turbulence-independent secondary "2" 
C     variables from the primary "2" variables.
C
C===================================================================70
      SUBROUTINE BLKIN

      use blpar_inc
      use xbl_inc
      IMPLICIT REAL(M)
      HVRAT = 0.0
C
C---- set edge Mach number ** 2
      M2    = U2*U2*HSTINV / (GM1BL*(1.0 - 0.5*U2*U2*HSTINV))
      TR2   = 1.0 + 0.5*GM1BL*M2
      M2_U2 = 2.0*M2*TR2/U2
      M2_MS = U2*U2*TR2    / (GM1BL*(1.0 - 0.5*U2*U2*HSTINV))
     &      * HSTINV_MS
C
C---- set edge static density (isentropic relation)
      R2    = RSTBL   *TR2**(-1.0/GM1BL)
      R2_U2 = -R2/TR2 * 0.5*M2_U2
      R2_MS = -R2/TR2 * 0.5*M2_MS
     &      + RSTBL_MS*TR2**(-1.0/GM1BL)
C
C---- set shape parameter
      H2    =  D2/T2
      H2_D2 = 1.0/T2
      H2_T2 = -H2/T2
C
C---- set edge static/stagnation enthalpy
      HERAT = 1.0 - 0.5*U2*U2*HSTINV
      HE_U2 =     -        U2*HSTINV
      HE_MS =     - 0.5*U2*U2*HSTINV_MS
C
C---- set molecular viscosity
      V2 = SQRT((HERAT)**3) * (1.0+HVRAT)/(HERAT+HVRAT)/REYBL
      V2_HE = V2*(1.5/HERAT - 1.0/(HERAT+HVRAT))
C
      V2_U2 =                        V2_HE*HE_U2
      V2_MS = -V2/REYBL * REYBL_MS + V2_HE*HE_MS
      V2_RE = -V2/REYBL * REYBL_RE
C
C---- set kinematic shape parameter
      CALL HKIN( H2, M2, HK2, HK2_H2, HK2_M2 )
C
      HK2_U2 =                HK2_M2*M2_U2
      HK2_T2 = HK2_H2*H2_T2
      HK2_D2 = HK2_H2*H2_D2
      HK2_MS =                HK2_M2*M2_MS
C
C---- set momentum thickness Reynolds number
      RT2    = R2*U2*T2/V2
      RT2_U2 = RT2*(1.0/U2 + R2_U2/R2 - V2_U2/V2)
      RT2_T2 = RT2/T2
      RT2_MS = RT2*(         R2_MS/R2 - V2_MS/V2)
      RT2_RE = RT2*(                  - V2_RE/V2)
C
      RETURN
      END ! BLKIN

C===================================================================70
C
C     Amplification rate routine for envelope e^n method.
C     Reference:
C                Drela, M., Giles, M.,
C               "Viscous/Inviscid Analysis of Transonic and
C                Low Reynolds Number Airfoils",
C                AIAA Journal, Oct. 1987.
C
C     NEW VERSION.   March 1991       (latest bug fix  July 93)
C          - m(H) correlation made more accurate up to H=20
C          - for H > 5, non-similar profiles are used 
C            instead of Falkner-Skan profiles.  These 
C            non-similar profiles have smaller reverse 
C            velocities, are more representative of typical 
C            separation bubble profiles.
C--------------------------------------------------------------
C
C     input :   HK     kinematic shape parameter
C               TH     momentum thickness
C               RT     momentum-thickness Reynolds number
C
C     output:   AX     envelope spatial amplification rate
C               AX_(.) sensitivity of AX to parameter (.)
C
C
C     Usage: The log of the envelope amplitude N(x) is
C            calculated by integrating AX (= dN/dx) with
C            respect to the streamwise distance x.
C                      x
C                     /
C              N(x) = | AX(H(x),Th(x),Rth(x)) dx
C                     /
C                      0
C            The integration can be started from the leading
C            edge since AX will be returned as zero when RT
C            is below the critical Rtheta.  Transition occurs
C            when N(x) reaches Ncrit (Ncrit= 9 is "standard").
C
C===================================================================70
      SUBROUTINE DAMPL( HK, TH, RT, AX, AX_HK, AX_TH, AX_RT )
      IMPLICIT REAL (A-H,M,O-Z)
ccc   DATA DGR / 0.04 /
      DATA DGR / 0.08 /
C
      HMI = 1.0/(HK - 1.0)
      HMI_HK = -HMI**2
C
C---- log10(Critical Rth) - H   correlation for Falkner-Skan profiles
      AA    = 2.492*HMI**0.43
      AA_HK =   (AA/HMI)*0.43 * HMI_HK
C
      BB    = TANH(14.0*HMI - 9.24)
      BB_HK = (1.0 - BB*BB) * 14.0 * HMI_HK
C
      GRCRIT = AA    + 0.7*(BB + 1.0)
      GRC_HK = AA_HK + 0.7* BB_HK
C
C
      GR = LOG10(RT)
      GR_RT = 1.0 / (2.3025851*RT)
C
      IF(GR .LT. GRCRIT-DGR) THEN
C
C----- no amplification for Rtheta < Rcrit
       AX    = 0.
       AX_HK = 0.
       AX_TH = 0.
       AX_RT = 0.
C
      ELSE
C
C----- Set steep cubic ramp used to turn on AX smoothly as Rtheta 
C-     exceeds Rcrit (previously, this was done discontinuously).
C-     The ramp goes between  -DGR < log10(Rtheta/Rcrit) < DGR
C
       RNORM = (GR - (GRCRIT-DGR)) / (2.0*DGR)
       RN_HK =     -  GRC_HK       / (2.0*DGR)
       RN_RT =  GR_RT              / (2.0*DGR)
C
       IF(RNORM .GE. 1.0) THEN
        RFAC    = 1.0
        RFAC_HK = 0.
        RFAC_RT = 0.
       ELSE
        RFAC    = 3.0*RNORM**2 - 2.0*RNORM**3
        RFAC_RN = 6.0*RNORM    - 6.0*RNORM**2
C
        RFAC_HK = RFAC_RN*RN_HK
        RFAC_RT = RFAC_RN*RN_RT
       ENDIF
C
C----- Amplification envelope slope correlation for Falkner-Skan
       ARG    = 3.87*HMI    - 2.52
       ARG_HK = 3.87*HMI_HK
C
       EX    = EXP(-ARG**2)
       EX_HK = EX * (-2.0*ARG*ARG_HK)
C
       DADR    = 0.028*(HK-1.0) - 0.0345*EX
       DADR_HK = 0.028          - 0.0345*EX_HK
C
C----- new m(H) correlation    1 March 91
       AF = -0.05 + 2.7*HMI -  5.5*HMI**2 + 3.0*HMI**3
       AF_HMI =     2.7     - 11.0*HMI    + 9.0*HMI**2
       AF_HK = AF_HMI*HMI_HK
C
       AX    = (AF   *DADR/TH                ) * RFAC
       AX_HK = (AF_HK*DADR/TH + AF*DADR_HK/TH) * RFAC
     &       + (AF   *DADR/TH                ) * RFAC_HK
       AX_TH = -AX/TH
       AX_RT = (AF   *DADR/TH                ) * RFAC_RT
C
      ENDIF
C
      RETURN
      END ! DAMPL

C===================================================================70 
C
C     Amplification rate routine for modified envelope e^n method.
C     Reference: 
C                Drela, M., Giles, M.,
C               "Viscous/Inviscid Analysis of Transonic and 
C                Low Reynolds Number Airfoils", 
C                AIAA Journal, Oct. 1987.
C
C     NEWER VERSION.   Nov 1996
C          - Amplification rate changes to the Orr-Sommerfeld
C              maximum ai(H,Rt) function for H > 4 .
C          - This implicitly assumes that the frequency range
C              (around w = 0.09 Ue/theta) which experiences this 
C              maximum amplification rate contains the currently
C              most-amplified frequency.
C--------------------------------------------------------------
C
C     input :   HK     kinematic shape parameter
C               TH     momentum thickness
C               RT     momentum-thickness Reynolds number
C
C     output:   AX     envelope spatial amplification rate
C               AX_(.) sensitivity of AX to parameter (.)
C
C
C     Usage: The log of the envelope amplitude N(x) is 
C            calculated by integrating AX (= dN/dx) with 
C            respect to the streamwise distance x.
C                      x
C                     /
C              N(x) = | AX(H(x),Th(x),Rth(x)) dx
C                     /
C                      0
C            The integration can be started from the leading
C            edge since AX will be returned as zero when RT
C            is below the critical Rtheta.  Transition occurs
C            when N(x) reaches Ncrit (Ncrit= 9 is "standard").
C
C===================================================================70 
      SUBROUTINE DAMPL2( HK, TH, RT, AX, AX_HK, AX_TH, AX_RT )
      IMPLICIT REAL (A-H,M,O-Z)
      DATA DGR / 0.08 /
      DATA HK1, HK2 / 3.5, 4.0 /
C
      HMI = 1.0/(HK - 1.0)
      HMI_HK = -HMI**2
C
C---- log10(Critical Rth) -- H   correlation for Falkner-Skan profiles
      AA    = 2.492*HMI**0.43
      AA_HK =   (AA/HMI)*0.43 * HMI_HK
C
      BB    = TANH(14.0*HMI - 9.24)
      BB_HK = (1.0 - BB*BB) * 14.0 * HMI_HK
C
      GRC = AA    + 0.7*(BB + 1.0)
      GRC_HK = AA_HK + 0.7* BB_HK
C
C
      GR = LOG10(RT)
      GR_RT = 1.0 / (2.3025851*RT)
C
      IF(GR .LT. GRC-DGR) THEN
C
C----- no amplification for Rtheta < Rcrit
       AX    = 0.
       AX_HK = 0.
       AX_TH = 0.
       AX_RT = 0.
C
      ELSE
C
C----- Set steep cubic ramp used to turn on AX smoothly as Rtheta 
C-     exceeds Rcrit (previously, this was done discontinuously).
C-     The ramp goes between  -DGR < log10(Rtheta/Rcrit) < DGR
C
       RNORM = (GR - (GRC-DGR)) / (2.0*DGR)
       RN_HK =     -  GRC_HK       / (2.0*DGR)
       RN_RT =  GR_RT              / (2.0*DGR)
C
       IF(RNORM .GE. 1.0) THEN
        RFAC    = 1.0
        RFAC_HK = 0.
        RFAC_RT = 0.
       ELSE
        RFAC    = 3.0*RNORM**2 - 2.0*RNORM**3
        RFAC_RN = 6.0*RNORM    - 6.0*RNORM**2
C
        RFAC_HK = RFAC_RN*RN_HK
        RFAC_RT = RFAC_RN*RN_RT
       ENDIF
C
C
C----- set envelope amplification rate with respect to Rtheta
C-       DADR = d(N)/d(Rtheta) = f(H)
C
       ARG    = 3.87*HMI    - 2.52
       ARG_HK = 3.87*HMI_HK
C
       EX    = EXP(-ARG**2)
       EX_HK = EX * (-2.0*ARG*ARG_HK)
C
       DADR    = 0.028*(HK-1.0) - 0.0345*EX
       DADR_HK = 0.028          - 0.0345*EX_HK
C
C
C----- set conversion factor from d/d(Rtheta) to d/dx
C-       AF = Theta d(Rtheta)/dx = f(H)
C
       BRG = -20.0*HMI
       AF = -0.05 + 2.7*HMI -  5.5*HMI**2 + 3.0*HMI**3 + 0.1*EXP(BRG)
       AF_HMI =     2.7     - 11.0*HMI    + 9.0*HMI**2 - 2.0*EXP(BRG)
       AF_HK = AF_HMI*HMI_HK
C
C
C----- set amplification rate with respect to x, 
C-     with RFAC shutting off amplification when below Rcrit
C
       AX    = (AF   *DADR/TH                ) * RFAC
       AX_HK = (AF_HK*DADR/TH + AF*DADR_HK/TH) * RFAC
     &       + (AF   *DADR/TH                ) * RFAC_HK
       AX_TH = -AX/TH
       AX_RT = (AF   *DADR/TH                ) * RFAC_RT
C
      ENDIF
C
      IF(HK .LT. HK1) RETURN
C
C---- non-envelope max-amplification correction for separated profiles
C
      HNORM = (HK - HK1) / (HK2 - HK1)
      HN_HK =       1.0  / (HK2 - HK1)
C
C---- set blending fraction HFAC = 0..1 over HK1 < HK < HK2
      IF(HNORM .GE. 1.0) THEN
       HFAC = 1.0
       HF_HK = 0.
      ELSE
       HFAC  =  3.0*HNORM**2 - 2.0*HNORM**3
       HF_HK = (6.0*HNORM    - 6.0*HNORM**2)*HN_HK
      ENDIF
C
C---- "normal" envelope amplification rate AX1
      AX1    = AX
      AX1_HK = AX_HK
      AX1_TH = AX_TH
      AX1_RT = AX_RT
C
C---- set modified amplification rate AX2
      GR0 = 0.30 + 0.35 * EXP(-0.15*(HK-5.0))
      GR0_HK =   - 0.35 * EXP(-0.15*(HK-5.0)) * 0.15
C
      TNR = TANH(1.2*(GR - GR0))
      TNR_RT =  (1.0 - TNR**2)*1.2*GR_RT
      TNR_HK = -(1.0 - TNR**2)*1.2*GR0_HK
C
      AX2    = (0.086*TNR    -     0.25/(HK-1.0)**1.5) / TH
      AX2_HK = (0.086*TNR_HK + 1.5*0.25/(HK-1.0)**2.5) / TH
      AX2_RT = (0.086*TNR_RT                         ) / TH
      AX2_TH = -AX2/TH
C
      IF(AX2 .LT. 0.0) THEN
       AX2    = 0.0
       AX2_HK = 0.
       AX2_RT = 0.
       AX2_TH = 0.
      ENDIF
C
C---- blend the two amplification rates
      AX    = HFAC*AX2    + (1.0 - HFAC)*AX1
      AX_HK = HFAC*AX2_HK + (1.0 - HFAC)*AX1_HK + HF_HK*(AX2-AX1)
      AX_RT = HFAC*AX2_RT + (1.0 - HFAC)*AX1_RT
      AX_TH = HFAC*AX2_TH + (1.0 - HFAC)*AX1_TH
C
      RETURN
      END ! DAMPL2

C===================================================================70
C
C     Returns average amplification AX over interval 1..2
C
C===================================================================70
      SUBROUTINE AXSET( HK1,    T1,    RT1,    A1,
     &                  HK2,    T2,    RT2,    A2,  ACRIT, IDAMPV,
     &           AX, AX_HK1, AX_T1, AX_RT1, AX_A1,
     &               AX_HK2, AX_T2, AX_RT2, AX_A2 )
C
cC==========================
cC---- 1st-order -- based on "1" quantities only
c      CALL DAMPL( HK1, T1, RT1, AX1, AX1_HK1, AX1_T1, AX1_RT1 )
c      AX2_HK2 = 0.0
c      AX2_T2  = 0.0
c      AX2_RT2 = 0.0
cC
c      AX1_A1 = 0.0
c      AX2_A2 = 0.0
cC
c      AX     = AX1
c      AX_AX1 = 1.0
c      AX_AX2 = 0.0
cC
c      ARG = MIN( 20.0*(ACRIT-A1) , 20.0 )
c      EXN    = EXP(-ARG)
c      EXN_A1 = 20.0*EXN
c      EXN_A2 = 0.
cC
c      DAX    = EXN   * 0.0004/T1
c      DAX_A1 = EXN_A1* 0.0004/T1
c      DAX_A2 = 0.
c      DAX_T1 = -DAX/T1
c      DAX_T2 = 0.
C
C==========================
C---- 2nd-order
      IF(IDAMPV.EQ.0) THEN
       CALL DAMPL( HK1, T1, RT1, AX1, AX1_HK1, AX1_T1, AX1_RT1 )
       CALL DAMPL( HK2, T2, RT2, AX2, AX2_HK2, AX2_T2, AX2_RT2 )
      ELSE
       CALL DAMPL2( HK1, T1, RT1, AX1, AX1_HK1, AX1_T1, AX1_RT1 )
       CALL DAMPL2( HK2, T2, RT2, AX2, AX2_HK2, AX2_T2, AX2_RT2 )
      ENDIF
C
CC---- simple-average version
C      AXA = 0.5*(AX1 + AX2)
C      IF(AXA .LE. 0.0) THEN
C       AXA = 0.0
C       AXA_AX1 = 0.0
C       AXA_AX2 = 0.0
C      ELSE
C       AXA_AX1 = 0.5
C       AXA_AX2 = 0.5
C      ENDIF
C
C---- rms-average version (seems a little better on coarse grids)
      AXSQ = 0.5*(AX1**2 + AX2**2)
      IF(AXSQ .LE. 0.0) THEN
       AXA = 0.0
       AXA_AX1 = 0.0
       AXA_AX2 = 0.0
      ELSE
       AXA = SQRT(AXSQ)
       AXA_AX1 = 0.5*AX1/AXA
       AXA_AX2 = 0.5*AX2/AXA
      ENDIF
C
C----- small additional term to ensure  dN/dx > 0  near  N = Ncrit
       ARG = MIN( 20.0*(ACRIT-0.5*(A1+A2)) , 20.0 )
       IF(ARG.LE.0.0) THEN
        EXN    = 1.0
CC      EXN_AC = 0.
        EXN_A1 = 0.
        EXN_A2 = 0.
       ELSE
        EXN    = EXP(-ARG)
CC      EXN_AC = -20.0    *EXN
        EXN_A1 =  20.0*0.5*EXN
        EXN_A2 =  20.0*0.5*EXN
       ENDIF
C
       DAX    = EXN    * 0.002/(T1+T2)
CC     DAX_AC = EXN_AC * 0.002/(T1+T2)
       DAX_A1 = EXN_A1 * 0.002/(T1+T2)
       DAX_A2 = EXN_A2 * 0.002/(T1+T2)
       DAX_T1 = -DAX/(T1+T2)
       DAX_T2 = -DAX/(T1+T2)
C
c
c        DAX    = 0.
c        DAX_A1 = 0.
c        DAX_A2 = 0.
c        DAX_AC = 0.
c        DAX_T1 = 0.
c        DAX_T2 = 0.
C==========================
C
      AX     = AXA             + DAX
C
      AX_HK1 = AXA_AX1*AX1_HK1
      AX_T1  = AXA_AX1*AX1_T1  + DAX_T1
      AX_RT1 = AXA_AX1*AX1_RT1
      AX_A1  =                   DAX_A1
C
      AX_HK2 = AXA_AX2*AX2_HK2
      AX_T2  = AXA_AX2*AX2_T2  + DAX_T2
      AX_RT2 = AXA_AX2*AX2_RT2
      AX_A2  =                   DAX_A2
C
      RETURN
      END

C===================================================================70
      SUBROUTINE TRCHEK2
C
C     New second-order version:  December 1994.
C
C     Checks if transition occurs in the current interval X1..X2.
C     If transition occurs, then set transition location XT, and 
C     its sensitivities to "1" and "2" variables.  If no transition, 
C     set amplification AMPL2.
C
C
C     Solves the implicit amplification equation for N2:
C
C       N2 - N1     N'(XT,NT) + N'(X1,N1)
C       -------  =  ---------------------
C       X2 - X1               2
C
C     In effect, a 2-point central difference is used between
C     X1..X2 (no transition), or X1..XT (transition).  The switch
C     is done by defining XT,NT in the equation above depending
C     on whether N2 exceeds Ncrit.
C
C  If N2<Ncrit:  NT=N2    , XT=X2                  (no transition)
C
C  If N2>Ncrit:  NT=Ncrit , XT=(Ncrit-N1)/(N2-N1)  (transition)
C
C
C
      use xfoil_inc, only : SILENT_MODE, XFOIL_FAIL
      use blpar_inc
      use xbl_inc
      DATA DAEPS / 5.0E-5 /
CCC   DATA DAEPS / 1.0D-12 /
C
C---- save variables and sensitivities at IBL ("2") for future restoration
C     DP mod: to remove need for EQUIVALENCE and COM2
      call store_c2sav()
C      DO 5 ICOM=1, NCOM
C        C2SAV(ICOM) = COM2(ICOM)
C    5 CONTINUE
C
C---- calculate average amplification rate AX over X1..X2 interval
      CALL AXSET( HK1,    T1,    RT1, AMPL1,
     &            HK2,    T2,    RT2, AMPL2,  AMCRIT, IDAMPV,
     &     AX, AX_HK1, AX_T1, AX_RT1, AX_A1,
     &         AX_HK2, AX_T2, AX_RT2, AX_A2 )
C
C---- set initial guess for iterate N2 (AMPL2) at X2
      AMPL2 = AMPL1 + AX*(X2-X1)
C
C---- solve implicit system for amplification AMPL2
      DO 100 ITAM=1, 30
C
C---- define weighting factors WF1,WF2 for defining "T" quantities from 1,2
C
      IF(AMPL2 .LE. AMCRIT) THEN
C------ there is no transition yet,  "T" is the same as "2"
        AMPLT    = AMPL2
        AMPLT_A2 = 1.0
        SFA    = 1.0
        SFA_A1 = 0.
        SFA_A2 = 0.
      ELSE
C------ there is transition in X1..X2, "T" is set from N1, N2
        AMPLT    = AMCRIT
        AMPLT_A2 = 0.
        SFA    = (AMPLT - AMPL1)/(AMPL2-AMPL1)
        SFA_A1 = ( SFA  - 1.0  )/(AMPL2-AMPL1)
        SFA_A2 = (      - SFA  )/(AMPL2-AMPL1)
      ENDIF
C
      IF(XIFORC.LT.X2) THEN
        SFX    = (XIFORC - X1 )/(X2-X1)
        SFX_X1 = (SFX    - 1.0)/(X2-X1)
        SFX_X2 = (       - SFX)/(X2-X1)
        SFX_XF =  1.0          /(X2-X1)
      ELSE
        SFX    = 1.0
        SFX_X1 = 0.
        SFX_X2 = 0.
        SFX_XF = 0.
      ENDIF
C
C---- set weighting factor from free or forced transition
      IF(SFA.LT.SFX) THEN
        WF2    = SFA
        WF2_A1 = SFA_A1
        WF2_A2 = SFA_A2
        WF2_X1 = 0.
        WF2_X2 = 0.
        WF2_XF = 0.
      ELSE
        WF2    = SFX
        WF2_A1 = 0.
        WF2_A2 = 0.
        WF2_X1 = SFX_X1
        WF2_X2 = SFX_X2
        WF2_XF = SFX_XF
      ENDIF
C
C
C=====================
CC---- 1st-order (based on "1" quantites only, for testing)
C      WF2    = 0.0
C      WF2_A1 = 0.0
C      WF2_A2 = 0.0
C      WF2_X1 = 0.0
C      WF2_X2 = 0.0
C      WF2_XF = 0.0
C=====================
C
      WF1    = 1.0 - WF2
      WF1_A1 =     - WF2_A1
      WF1_A2 =     - WF2_A2
      WF1_X1 =     - WF2_X1
      WF1_X2 =     - WF2_X2
      WF1_XF =     - WF2_XF
C
C---- interpolate BL variables to XT
      XT    = X1*WF1    + X2*WF2
      TT    = T1*WF1    + T2*WF2
      DT    = D1*WF1    + D2*WF2
      UT    = U1*WF1    + U2*WF2
C
      XT_A2 = X1*WF1_A2 + X2*WF2_A2
      TT_A2 = T1*WF1_A2 + T2*WF2_A2
      DT_A2 = D1*WF1_A2 + D2*WF2_A2
      UT_A2 = U1*WF1_A2 + U2*WF2_A2
C
C---- temporarily set "2" variables from "T" for BLKIN
      X2 = XT
      T2 = TT
      D2 = DT
      U2 = UT
C
C---- calculate laminar secondary "T" variables HKT, RTT
      CALL BLKIN
C
      HKT    = HK2
      HKT_TT = HK2_T2
      HKT_DT = HK2_D2
      HKT_UT = HK2_U2
      HKT_MS = HK2_MS
C
      RTT    = RT2
      RTT_TT = RT2_T2
      RTT_UT = RT2_U2
      RTT_MS = RT2_MS
      RTT_RE = RT2_RE
C
C---- restore clobbered "2" variables, except for AMPL2
      AMSAVE = AMPL2
C     DP mod: to remove need for EQUIVALENCE and COM2
      call from_c2sav()
C      DO 8 ICOM=1, NCOM
C        COM2(ICOM) = C2SAV(ICOM)
C 8    CONTINUE
      AMPL2 = AMSAVE
C
C---- calculate amplification rate AX over current X1-XT interval
      CALL AXSET( HK1,    T1,    RT1, AMPL1,
     &            HKT,    TT,    RTT, AMPLT,  AMCRIT, IDAMPV,
     &     AX, AX_HK1, AX_T1, AX_RT1, AX_A1,
     &         AX_HKT, AX_TT, AX_RTT, AX_AT )
C
C---- punch out early if there is no amplification here
      IF(AX .LE. 0.0) GO TO 101
C
C---- set sensitivity of AX(A2)
      AX_A2 = (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_A2
     &      + (AX_HKT*HKT_DT                        )*DT_A2
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_A2
     &      +  AX_AT                                 *AMPLT_A2
C
C---- residual for implicit AMPL2 definition (amplification equation)
      RES    = AMPL2 - AMPL1 - AX   *(X2-X1) 
      RES_A2 = 1.0           - AX_A2*(X2-X1)
C
      DA2 = -RES/RES_A2
C
      RLX = 1.0
      DXT = XT_A2*DA2
C
      IF(RLX*ABS(DXT/(X2-X1)) .GT. 0.05) RLX = 0.05*ABS((X2-X1)/DXT)
      IF(RLX*ABS(DA2)         .GT. 1.0 ) RLX = 1.0 *ABS(   1.0 /DA2)
C
C---- check if converged
      IF(ABS(DA2) .LT. DAEPS) GO TO 101
C
      IF((AMPL2.GT.AMCRIT .AND. AMPL2+RLX*DA2.LT.AMCRIT).OR.
     &   (AMPL2.LT.AMCRIT .AND. AMPL2+RLX*DA2.GT.AMCRIT)    ) THEN
C------ limited Newton step so AMPL2 doesn't step across AMCRIT either way
        AMPL2 = AMCRIT
      ELSE
C------ regular Newton step
        AMPL2 = AMPL2 + RLX*DA2
      ENDIF
C
 100  CONTINUE
C     DP mod: check for infinite loop condition
      IF (ISNAN(AMPL1) .OR. ISNAN(AMPLT) .OR. ISNAN(AMPL2)) THEN
        XFOIL_FAIL = .TRUE.
        RETURN
      ENDIF

C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE) THEN
        WRITE(*,*) 'TRCHEK2: N2 convergence failed.'
        WRITE(*,6700) X1, XT, X2, AMPL1, AMPLT, AMPL2, AX, DA2
      ENDIF
 6700 FORMAT(1X,'x:', 3F9.5,'  N:',3F7.3,'  Nx:',F8.3,'   dN:',E10.3)
C
 101  CONTINUE
C
C
C---- test for free or forced transition
      TRFREE = AMPL2 .GE. AMCRIT
      TRFORC = XIFORC.GT.X1 .AND. XIFORC.LE.X2
C
C---- set transition interval flag
      TRAN = TRFORC .OR. TRFREE
C
      IF(.NOT.TRAN) RETURN
C
C---- resolve if both forced and free transition
      IF(TRFREE .AND. TRFORC) THEN
       TRFORC = XIFORC .LT. XT
       TRFREE = XIFORC .GE. XT
      ENDIF
C
      IF(TRFORC) THEN
C----- if forced transition, then XT is prescribed,
C-     no sense calculating the sensitivities, since we know them...
       XT = XIFORC
       XT_A1 = 0.
       XT_X1 = 0.
       XT_T1 = 0.
       XT_D1 = 0.
       XT_U1 = 0.
       XT_X2 = 0.
       XT_T2 = 0.
       XT_D2 = 0.
       XT_U2 = 0.
       XT_MS = 0.
       XT_RE = 0.
       XT_XF = 1.0
       RETURN
      ENDIF
C
C---- free transition ... set sensitivities of XT
C
C---- XT( X1 X2 A1 A2 XF ),  TT( T1 T2 A1 A2 X1 X2 XF),   DT( ...
CC    XT    = X1*WF1    + X2*WF2
CC    TT    = T1*WF1    + T2*WF2
CC    DT    = D1*WF1    + D2*WF2
CC    UT    = U1*WF1    + U2*WF2
C
      XT_X1 =    WF1
      TT_T1 =    WF1
      DT_D1 =    WF1
      UT_U1 =    WF1
C
      XT_X2 =                WF2
      TT_T2 =                WF2
      DT_D2 =                WF2
      UT_U2 =                WF2
C
      XT_A1 = X1*WF1_A1 + X2*WF2_A1
      TT_A1 = T1*WF1_A1 + T2*WF2_A1
      DT_A1 = D1*WF1_A1 + D2*WF2_A1
      UT_A1 = U1*WF1_A1 + U2*WF2_A1
C
CC    XT_A2 = X1*WF1_A2 + X2*WF2_A2
CC    TT_A2 = T1*WF1_A2 + T2*WF2_A2
CC    DT_A2 = D1*WF1_A2 + D2*WF2_A2
CC    UT_A2 = U1*WF1_A2 + U2*WF2_A2
C
      XT_X1 = X1*WF1_X1 + X2*WF2_X1 + XT_X1
      TT_X1 = T1*WF1_X1 + T2*WF2_X1
      DT_X1 = D1*WF1_X1 + D2*WF2_X1
      UT_X1 = U1*WF1_X1 + U2*WF2_X1
C
      XT_X2 = X1*WF1_X2 + X2*WF2_X2 + XT_X2
      TT_X2 = T1*WF1_X2 + T2*WF2_X2
      DT_X2 = D1*WF1_X2 + D2*WF2_X2
      UT_X2 = U1*WF1_X2 + U2*WF2_X2
C
      XT_XF = X1*WF1_XF + X2*WF2_XF
      TT_XF = T1*WF1_XF + T2*WF2_XF
      DT_XF = D1*WF1_XF + D2*WF2_XF
      UT_XF = U1*WF1_XF + U2*WF2_XF
C
C---- at this point, AX = AX( HK1, T1, RT1, A1, HKT, TT, RTT, AT )
C
C---- set sensitivities of AX( T1 D1 U1 A1 T2 D2 U2 A2 MS RE )
      AX_T1 =  AX_HK1*HK1_T1 + AX_T1 + AX_RT1*RT1_T1
     &      + (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_T1
      AX_D1 =  AX_HK1*HK1_D1
     &      + (AX_HKT*HKT_DT                        )*DT_D1
      AX_U1 =  AX_HK1*HK1_U1         + AX_RT1*RT1_U1
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_U1
      AX_A1 =  AX_A1
     &      + (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_A1
     &      + (AX_HKT*HKT_DT                        )*DT_A1
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_A1
      AX_X1 = (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_X1
     &      + (AX_HKT*HKT_DT                        )*DT_X1
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_X1
C
      AX_T2 = (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_T2
      AX_D2 = (AX_HKT*HKT_DT                        )*DT_D2
      AX_U2 = (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_U2
      AX_A2 =  AX_AT                                 *AMPLT_A2
     &      + (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_A2
     &      + (AX_HKT*HKT_DT                        )*DT_A2
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_A2
      AX_X2 = (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_X2
     &      + (AX_HKT*HKT_DT                        )*DT_X2
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_X2
C
      AX_XF = (AX_HKT*HKT_TT + AX_TT + AX_RTT*RTT_TT)*TT_XF
     &      + (AX_HKT*HKT_DT                        )*DT_XF
     &      + (AX_HKT*HKT_UT         + AX_RTT*RTT_UT)*UT_XF
C
      AX_MS =  AX_HKT*HKT_MS         + AX_RTT*RTT_MS
     &      +  AX_HK1*HK1_MS         + AX_RT1*RT1_MS
      AX_RE =                          AX_RTT*RTT_RE
     &                               + AX_RT1*RT1_RE
C
C
C---- set sensitivities of residual RES
CCC   RES  = AMPL2 - AMPL1 - AX*(X2-X1)
      Z_AX =               -    (X2-X1)
C
      Z_A1 = Z_AX*AX_A1 - 1.0
      Z_T1 = Z_AX*AX_T1
      Z_D1 = Z_AX*AX_D1
      Z_U1 = Z_AX*AX_U1
      Z_X1 = Z_AX*AX_X1 + AX
C
      Z_A2 = Z_AX*AX_A2 + 1.0
      Z_T2 = Z_AX*AX_T2
      Z_D2 = Z_AX*AX_D2
      Z_U2 = Z_AX*AX_U2
      Z_X2 = Z_AX*AX_X2 - AX
C
      Z_XF = Z_AX*AX_XF
      Z_MS = Z_AX*AX_MS
      Z_RE = Z_AX*AX_RE
C
C---- set sensitivities of XT, with RES being stationary for A2 constraint
      XT_A1 = XT_A1 - (XT_A2/Z_A2)*Z_A1
      XT_T1 =       - (XT_A2/Z_A2)*Z_T1
      XT_D1 =       - (XT_A2/Z_A2)*Z_D1
      XT_U1 =       - (XT_A2/Z_A2)*Z_U1
      XT_X1 = XT_X1 - (XT_A2/Z_A2)*Z_X1
      XT_T2 =       - (XT_A2/Z_A2)*Z_T2
      XT_D2 =       - (XT_A2/Z_A2)*Z_D2
      XT_U2 =       - (XT_A2/Z_A2)*Z_U2
      XT_X2 = XT_X2 - (XT_A2/Z_A2)*Z_X2
      XT_MS =       - (XT_A2/Z_A2)*Z_MS
      XT_RE =       - (XT_A2/Z_A2)*Z_RE
      XT_XF = 0.0
C
      RETURN
      END

C===================================================================70
C 1st-order amplification equation
C===================================================================70
      SUBROUTINE TRCHEK
C
cc      CALL TRCHEK1
C
C---- 2nd-order amplification equation
      CALL TRCHEK2
C
      RETURN
      END

C===================================================================70
C===================================================================70
      SUBROUTINE HCT( HK, MSQ, HC, HC_HK, HC_MSQ )
      REAL MSQ
C
C---- density shape parameter    (from Whitfield)
      HC     = MSQ * (0.064/(HK-0.8) + 0.251)
      HC_HK  = MSQ * (-.064/(HK-0.8)**2     )
      HC_MSQ =        0.064/(HK-0.8) + 0.251
C
      RETURN
      END

C===================================================================70
C===================================================================70
      SUBROUTINE HSL( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
      REAL MSQ
C
C---- Laminar HS correlation
      IF(HK.LT.4.35) THEN
       TMP = HK - 4.35
       HS    = 0.0111*TMP**2/(HK+1.0)
     &       - 0.0278*TMP**3/(HK+1.0)  + 1.528
     &       - 0.0002*(TMP*HK)**2
       HS_HK = 0.0111*(2.0*TMP    - TMP**2/(HK+1.0))/(HK+1.0)
     &       - 0.0278*(3.0*TMP**2 - TMP**3/(HK+1.0))/(HK+1.0)
     &       - 0.0002*2.0*TMP*HK * (TMP + HK)
      ELSE
       HS    = 0.015*    (HK-4.35)**2/HK + 1.528
       HS_HK = 0.015*2.0*(HK-4.35)   /HK
     &       - 0.015*    (HK-4.35)**2/HK**2
      ENDIF
C
      HS_RT  = 0.
      HS_MSQ = 0.
C
      RETURN
      END

C===================================================================70
C---- Turbulent HS correlation
C===================================================================70
      SUBROUTINE HST( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
      IMPLICIT REAL (A-H,M,O-Z)
C
C
      DATA HSMIN, DHSINF / 1.500, 0.015 /
C
C---- ###  12/4/94
C---- limited Rtheta dependence for Rtheta < 200
C
C
      IF(RT.GT.400.0) THEN
       HO    = 3.0 + 400.0/RT
       HO_RT =     - 400.0/RT**2
      ELSE
       HO    = 4.0
       HO_RT = 0.
      ENDIF
C
      IF(RT.GT.200.0) THEN
       RTZ    = RT
       RTZ_RT = 1.
      ELSE
       RTZ    = 200.0
       RTZ_RT = 0.
      ENDIF
C
      IF(HK.LT.HO) THEN
C----- attached branch
C=======================================================
C----- old correlation
C-     (from Swafford profiles)
c       SRT = SQRT(RT)
c       HEX = (HO-HK)**1.6
c       RTMP = 0.165 - 1.6/SRT
c       HS    = HSMIN + 4.0/RT + RTMP*HEX/HK
c       HS_HK = RTMP*HEX/HK*(-1.6/(HO-HK) - 1.0/HK)
c       HS_RT = -4.0/RT**2 + HEX/HK*0.8/SRT/RT
c     &             + RTMP*HEX/HK*1.6/(HO-HK)*HO_RT
C=======================================================
C----- new correlation  29 Nov 91
C-     (from  arctan(y+) + Schlichting  profiles)
       HR    = ( HO - HK)/(HO-1.0)
       HR_HK =      - 1.0/(HO-1.0)
       HR_RT = (1.0 - HR)/(HO-1.0) * HO_RT
       HS    = (2.0-HSMIN-4.0/RTZ)*HR**2  * 1.5/(HK+0.5) + HSMIN
     &       + 4.0/RTZ
       HS_HK =-(2.0-HSMIN-4.0/RTZ)*HR**2  * 1.5/(HK+0.5)**2
     &       + (2.0-HSMIN-4.0/RTZ)*HR*2.0 * 1.5/(HK+0.5) * HR_HK
       HS_RT = (2.0-HSMIN-4.0/RTZ)*HR*2.0 * 1.5/(HK+0.5) * HR_RT
     &       + (HR**2 * 1.5/(HK+0.5) - 1.0)*4.0/RTZ**2 * RTZ_RT
C
      ELSE
C
C----- separated branch
       GRT = LOG(RTZ)
       HDIF = HK - HO 
       RTMP = HK - HO + 4.0/GRT
       HTMP    = 0.007*GRT/RTMP**2 + DHSINF/HK
       HTMP_HK = -.014*GRT/RTMP**3 - DHSINF/HK**2
       HTMP_RT = -.014*GRT/RTMP**3 * (-HO_RT - 4.0/GRT**2/RTZ * RTZ_RT)
     &         + 0.007    /RTMP**2 / RTZ * RTZ_RT
       HS    = HDIF**2 * HTMP + HSMIN + 4.0/RTZ
       HS_HK = HDIF*2.0* HTMP
     &       + HDIF**2 * HTMP_HK
       HS_RT = HDIF**2 * HTMP_RT      - 4.0/RTZ**2 * RTZ_RT
     &       + HDIF*2.0* HTMP * (-HO_RT)
C
      ENDIF
C
C---- fudge HS slightly to make sure   HS -> 2   as   HK -> 1
C-    (unnecessary with new correlation)
c      HTF    = 0.485/9.0 * (HK-4.0)**2/HK  +  1.515
c      HTF_HK = 0.485/9.0 * (1.0-16.0/HK**2)
c      ARG = MAX( 10.0*(1.0 - HK) , -15.0 )
c      HXX = EXP(ARG)
c      HXX_HK = -10.0*HXX
cC
c      HS_HK  = (1.0-HXX)*HS_HK  +  HXX*HTF_HK
c     &       + (        -HS     +      HTF    )*HXX_HK
c      HS_RT  = (1.0-HXX)*HS_RT
c      HS     = (1.0-HXX)*HS     +  HXX*HTF
C
C---- Whitfield's minor additional compressibility correction
      FM = 1.0 + 0.014*MSQ
      HS     = ( HS + 0.028*MSQ ) / FM
      HS_HK  = ( HS_HK          ) / FM
      HS_RT  = ( HS_RT          ) / FM
      HS_MSQ = 0.028/FM  -  0.014*HS/FM
C
      RETURN
      END

C===================================================================70
C---- Laminar skin friction function  ( Cf )    ( from Falkner-Skan )
C===================================================================70
      SUBROUTINE CFL( HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ )
      REAL*8 :: MSQ
C
      IF(HK.LT.5.5) THEN
       TMP = (5.5-HK)**3 / (HK+1.0)
       CF    = ( 0.0727*TMP                      - 0.07       )/RT
       CF_HK = ( -.0727*TMP*3.0/(5.5-HK) - 0.0727*TMP/(HK+1.0))/RT
      ELSE
       TMP = 1.0 - 1.0/(HK-4.5)
       CF    = ( 0.015*TMP**2      - 0.07  ) / RT
       CF_HK = ( 0.015*TMP*2.0/(HK-4.5)**2 ) / RT
      ENDIF
      CF_RT = -CF/RT
      CF_MSQ = 0.0
C
      RETURN
      END

C===================================================================70
C---- Turbulent skin friction function  ( Cf )    (Coles)
C===================================================================70
      SUBROUTINE CFT( HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ )

      use blpar_inc
      IMPLICIT REAL (A-H,M,O-Z)
C
      DATA GAM /1.4/
C
      GM1 = GAM - 1.0
      FC = SQRT(1.0 + 0.5*GM1*MSQ)
      GRT = LOG(RT/FC)
      GRT = MAX(GRT,3.0)
C
      GEX = -1.74 - 0.31*HK
C
      ARG = -1.33*HK
      ARG = MAX(-20.0, ARG )
C
      THK = TANH(4.0 - HK/0.875)
C
      CFO =  CFFAC * 0.3*EXP(ARG) * (GRT/2.3026)**GEX
      CF     = ( CFO  +  1.1E-4*(THK-1.0) ) / FC
      CF_HK  = (-1.33*CFO - 0.31*LOG(GRT/2.3026)*CFO
     &         - 1.1E-4*(1.0-THK**2) / 0.875    ) / FC
      CF_RT  = GEX*CFO/(FC*GRT) / RT
      CF_MSQ = GEX*CFO/(FC*GRT) * (-0.25*GM1/FC**2) - 0.25*GM1*CF/FC**2
C
      RETURN
      END ! CFT

C===================================================================70
C---- Laminar dissipation function  ( 2 CD/H* )     (from Falkner-Skan)
C===================================================================70
      SUBROUTINE DIL( HK, RT, DI, DI_HK, DI_RT )
C
      IF(HK.LT.4.0) THEN
       DI    = ( 0.00205  *  (4.0-HK)**5.5 + 0.207 ) / RT
       DI_HK = ( -.00205*5.5*(4.0-HK)**4.5         ) / RT
      ELSE
       HKB = HK - 4.0
       DEN = 1.0 + 0.02*HKB**2
       DI    = ( -.0016  *  HKB**2  /DEN   + 0.207             ) / RT
       DI_HK = ( -.0016*2.0*HKB*(1.0/DEN - 0.02*HKB**2/DEN**2) ) / RT
      ENDIF
      DI_RT = -DI/RT
C
      RETURN
      END

C===================================================================70
C===================================================================70
      SUBROUTINE DILW( HK, RT, DI, DI_HK, DI_RT )
      REAL MSQ
C
      MSQ = 0.
      CALL HSL( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
C
C---- Laminar wake dissipation function  ( 2 CD/H* )
      RCD    =  1.10 * (1.0 - 1.0/HK)**2  / HK
      RCD_HK = -1.10 * (1.0 - 1.0/HK)*2.0 / HK**3
     &       - RCD/HK
C
      DI    = 2.0*RCD   /(HS*RT)
      DI_HK = 2.0*RCD_HK/(HS*RT) - (DI/HS)*HS_HK
      DI_RT = -DI/RT             - (DI/HS)*HS_RT
C
      RETURN
      END

C===================================================================70
C
C     Calculates all secondary "2" variables from
C     the primary "2" variables X2, U2, T2, D2, S2.
C     Also calculates the sensitivities of the
C     secondary variables wrt the primary variables.
C
C      ITYP = 1 :  laminar
C      ITYP = 2 :  turbulent
C      ITYP = 3 :  turbulent wake
C
C===================================================================70
      SUBROUTINE BLVAR(ITYP)

      use blpar_inc
      use xbl_inc
      IMPLICIT REAL(M)
C
      IF(ITYP.EQ.3) HK2 = MAX(HK2,1.00005)
      IF(ITYP.NE.3) HK2 = MAX(HK2,1.05000)
C
C---- density thickness shape parameter     ( H** )
      CALL HCT( HK2, M2, HC2, HC2_HK2, HC2_M2 )
      HC2_U2 = HC2_HK2*HK2_U2 + HC2_M2*M2_U2
      HC2_T2 = HC2_HK2*HK2_T2
      HC2_D2 = HC2_HK2*HK2_D2
      HC2_MS = HC2_HK2*HK2_MS + HC2_M2*M2_MS
C
C---- set KE thickness shape parameter from  H - H*  correlations
      IF(ITYP.EQ.1) THEN
       CALL HSL( HK2, RT2, M2, HS2, HS2_HK2, HS2_RT2, HS2_M2 )
      ELSE
       CALL HST( HK2, RT2, M2, HS2, HS2_HK2, HS2_RT2, HS2_M2 )
      ENDIF
C
      HS2_U2 = HS2_HK2*HK2_U2 + HS2_RT2*RT2_U2 + HS2_M2*M2_U2
      HS2_T2 = HS2_HK2*HK2_T2 + HS2_RT2*RT2_T2
      HS2_D2 = HS2_HK2*HK2_D2
      HS2_MS = HS2_HK2*HK2_MS + HS2_RT2*RT2_MS + HS2_M2*M2_MS
      HS2_RE =                  HS2_RT2*RT2_RE
C
C---- normalized slip velocity  Us
      US2     = 0.5*HS2*( 1.0 - (HK2-1.0)/(GBCON*H2) )
      US2_HS2 = 0.5  *  ( 1.0 - (HK2-1.0)/(GBCON*H2) )
      US2_HK2 = 0.5*HS2*(     -  1.0     /(GBCON*H2) )
      US2_H2  = 0.5*HS2*        (HK2-1.0)/(GBCON*H2**2)
C
      US2_U2 = US2_HS2*HS2_U2 + US2_HK2*HK2_U2
      US2_T2 = US2_HS2*HS2_T2 + US2_HK2*HK2_T2 + US2_H2*H2_T2
      US2_D2 = US2_HS2*HS2_D2 + US2_HK2*HK2_D2 + US2_H2*H2_D2
      US2_MS = US2_HS2*HS2_MS + US2_HK2*HK2_MS
      US2_RE = US2_HS2*HS2_RE
C
      IF(ITYP.LE.2 .AND. US2.GT.0.95) THEN
CCC       WRITE(*,*) 'BLVAR: Us clamped:', US2
       US2 = 0.98
       US2_U2 = 0.
       US2_T2 = 0.
       US2_D2 = 0.
       US2_MS = 0.
       US2_RE = 0.
      ENDIF
C
      IF(ITYP.EQ.3 .AND. US2.GT.0.99995) THEN
CCC       WRITE(*,*) 'BLVAR: Wake Us clamped:', US2
       US2 = 0.99995
       US2_U2 = 0.
       US2_T2 = 0.
       US2_D2 = 0.
       US2_MS = 0.
       US2_RE = 0.
      ENDIF
C
C---- equilibrium wake layer shear coefficient (Ctau)EQ ** 1/2
C   ...  NEW  12 Oct 94
      GCC = 0.0
      HKC = HK2 - 1.0
      HKC_HK2 = 1.0
      HKC_RT2 = 0.0
      IF(ITYP.EQ.2) THEN
       GCC = GCCON
       HKC     = HK2 - 1.0 - GCC/RT2
       HKC_HK2 = 1.0
       HKC_RT2 =             GCC/RT2**2
       IF(HKC .LT. 0.01) THEN
        HKC = 0.01
        HKC_HK2 = 0.0
        HKC_RT2 = 0.0
       ENDIF
      ENDIF
C
      HKB = HK2 - 1.0
      USB = 1.0 - US2
      CQ2     =
     &    SQRT( CTCON*HS2*HKB*HKC**2 / (USB*H2*HK2**2) )
      CQ2_HS2 = CTCON    *HKB*HKC**2 / (USB*H2*HK2**2)       * 0.5/CQ2
      CQ2_US2 = CTCON*HS2*HKB*HKC**2 / (USB*H2*HK2**2) / USB * 0.5/CQ2
      CQ2_HK2 = CTCON*HS2    *HKC**2 / (USB*H2*HK2**2)       * 0.5/CQ2
     &        - CTCON*HS2*HKB*HKC**2 / (USB*H2*HK2**3) * 2.0 * 0.5/CQ2
     &        + CTCON*HS2*HKB*HKC    / (USB*H2*HK2**2) * 2.0 * 0.5/CQ2
     &         *HKC_HK2
      CQ2_RT2 = CTCON*HS2*HKB*HKC    / (USB*H2*HK2**2) * 2.0 * 0.5/CQ2
     &         *HKC_RT2
      CQ2_H2  =-CTCON*HS2*HKB*HKC**2 / (USB*H2*HK2**2) / H2  * 0.5/CQ2
C
      CQ2_U2 = CQ2_HS2*HS2_U2 + CQ2_US2*US2_U2 + CQ2_HK2*HK2_U2
      CQ2_T2 = CQ2_HS2*HS2_T2 + CQ2_US2*US2_T2 + CQ2_HK2*HK2_T2
      CQ2_D2 = CQ2_HS2*HS2_D2 + CQ2_US2*US2_D2 + CQ2_HK2*HK2_D2
      CQ2_MS = CQ2_HS2*HS2_MS + CQ2_US2*US2_MS + CQ2_HK2*HK2_MS
      CQ2_RE = CQ2_HS2*HS2_RE + CQ2_US2*US2_RE
C
      CQ2_U2 = CQ2_U2                + CQ2_RT2*RT2_U2
      CQ2_T2 = CQ2_T2 + CQ2_H2*H2_T2 + CQ2_RT2*RT2_T2
      CQ2_D2 = CQ2_D2 + CQ2_H2*H2_D2
      CQ2_MS = CQ2_MS                + CQ2_RT2*RT2_MS
      CQ2_RE = CQ2_RE                + CQ2_RT2*RT2_RE
C
C
C---- set skin friction coefficient 
      IF(ITYP.EQ.3) THEN
C----- wake
       CF2     = 0.
       CF2_HK2 = 0.
       CF2_RT2 = 0.
       CF2_M2  = 0.
      ELSE IF(ITYP.EQ.1) THEN
C----- laminar
       CALL CFL( HK2, RT2, M2, CF2, CF2_HK2, CF2_RT2, CF2_M2 )
      ELSE
C----- turbulent
       CALL CFT( HK2, RT2, M2, CF2, CF2_HK2, CF2_RT2, CF2_M2 )
       CALL CFL( HK2, RT2, M2, CF2L,CF2L_HK2,CF2L_RT2,CF2L_M2)
       IF(CF2L.GT.CF2) THEN
C------- laminar Cf is greater than turbulent Cf -- use laminar
C-       (this will only occur for unreasonably small Rtheta)
ccc      write(*,*) 'Cft Cfl Rt Hk:', CF2, CF2L, RT2, HK2, X2
         CF2     = CF2L
         CF2_HK2 = CF2L_HK2
         CF2_RT2 = CF2L_RT2
         CF2_M2  = CF2L_M2
       ENDIF
      ENDIF
C
      CF2_U2 = CF2_HK2*HK2_U2 + CF2_RT2*RT2_U2 + CF2_M2*M2_U2
      CF2_T2 = CF2_HK2*HK2_T2 + CF2_RT2*RT2_T2
      CF2_D2 = CF2_HK2*HK2_D2
      CF2_MS = CF2_HK2*HK2_MS + CF2_RT2*RT2_MS + CF2_M2*M2_MS
      CF2_RE =                  CF2_RT2*RT2_RE
C
C---- dissipation function    2 CD / H*
      IF(ITYP.EQ.1) THEN
C
C----- laminar
       CALL DIL( HK2, RT2, DI2, DI2_HK2, DI2_RT2 )
C
       DI2_U2 = DI2_HK2*HK2_U2 + DI2_RT2*RT2_U2
       DI2_T2 = DI2_HK2*HK2_T2 + DI2_RT2*RT2_T2
       DI2_D2 = DI2_HK2*HK2_D2
       DI2_S2 = 0.
       DI2_MS = DI2_HK2*HK2_MS + DI2_RT2*RT2_MS
       DI2_RE =                  DI2_RT2*RT2_RE
C
      ELSE IF(ITYP.EQ.2) THEN
C
CCC       CALL DIT(     HS2,     US2,     CF2,     S2, DI2,
CCC     &           DI2_HS2, DI2_US2, DI2_CF2, DI2_S2      )
C
C----- turbulent wall contribution
       CALL CFT(HK2, RT2, M2, CF2T, CF2T_HK2, CF2T_RT2, CF2T_M2)
       CF2T_U2 = CF2T_HK2*HK2_U2 + CF2T_RT2*RT2_U2 + CF2T_M2*M2_U2
       CF2T_T2 = CF2T_HK2*HK2_T2 + CF2T_RT2*RT2_T2
       CF2T_D2 = CF2T_HK2*HK2_D2
       CF2T_MS = CF2T_HK2*HK2_MS + CF2T_RT2*RT2_MS + CF2T_M2*M2_MS
       CF2T_RE =                   CF2T_RT2*RT2_RE
C
       DI2      =  ( 0.5*CF2T*US2 ) * 2.0/HS2
       DI2_HS2  = -( 0.5*CF2T*US2 ) * 2.0/HS2**2
       DI2_US2  =  ( 0.5*CF2T     ) * 2.0/HS2
       DI2_CF2T =  ( 0.5     *US2 ) * 2.0/HS2
C
       DI2_S2 = 0.0
       DI2_U2 = DI2_HS2*HS2_U2 + DI2_US2*US2_U2 + DI2_CF2T*CF2T_U2
       DI2_T2 = DI2_HS2*HS2_T2 + DI2_US2*US2_T2 + DI2_CF2T*CF2T_T2
       DI2_D2 = DI2_HS2*HS2_D2 + DI2_US2*US2_D2 + DI2_CF2T*CF2T_D2
       DI2_MS = DI2_HS2*HS2_MS + DI2_US2*US2_MS + DI2_CF2T*CF2T_MS
       DI2_RE = DI2_HS2*HS2_RE + DI2_US2*US2_RE + DI2_CF2T*CF2T_RE
C
C
C----- set minimum Hk for wake layer to still exist
       GRT = LOG(RT2)
       HMIN = 1.0 + 2.1/GRT
       HM_RT2 = -(2.1/GRT**2) / RT2
C
C----- set factor DFAC for correcting wall dissipation for very low Hk
       FL = (HK2-1.0)/(HMIN-1.0)
       FL_HK2 =   1.0/(HMIN-1.0)
       FL_RT2 = ( -FL/(HMIN-1.0) ) * HM_RT2
C
       TFL = TANH(FL)
       DFAC  = 0.5 + 0.5* TFL
       DF_FL =       0.5*(1.0 - TFL**2)
C
       DF_HK2 = DF_FL*FL_HK2
       DF_RT2 = DF_FL*FL_RT2
C
       DI2_S2 = DI2_S2*DFAC
       DI2_U2 = DI2_U2*DFAC + DI2*(DF_HK2*HK2_U2 + DF_RT2*RT2_U2)
       DI2_T2 = DI2_T2*DFAC + DI2*(DF_HK2*HK2_T2 + DF_RT2*RT2_T2)
       DI2_D2 = DI2_D2*DFAC + DI2*(DF_HK2*HK2_D2                )
       DI2_MS = DI2_MS*DFAC + DI2*(DF_HK2*HK2_MS + DF_RT2*RT2_MS)
       DI2_RE = DI2_RE*DFAC + DI2*(                DF_RT2*RT2_RE)
       DI2    = DI2   *DFAC
C
      ELSE
C
C----- zero wall contribution for wake
       DI2    = 0.0
       DI2_S2 = 0.0
       DI2_U2 = 0.0
       DI2_T2 = 0.0
       DI2_D2 = 0.0
       DI2_MS = 0.0
       DI2_RE = 0.0
C
      ENDIF
C
C
C---- Add on turbulent outer layer contribution
      IF(ITYP.NE.1) THEN
C
       DD     =  S2**2 * (0.995-US2) * 2.0/HS2
       DD_HS2 = -S2**2 * (0.995-US2) * 2.0/HS2**2
       DD_US2 = -S2**2               * 2.0/HS2
       DD_S2  =  S2*2.0* (0.995-US2) * 2.0/HS2
C
       DI2    = DI2    + DD
       DI2_S2 =          DD_S2
       DI2_U2 = DI2_U2 + DD_HS2*HS2_U2 + DD_US2*US2_U2
       DI2_T2 = DI2_T2 + DD_HS2*HS2_T2 + DD_US2*US2_T2
       DI2_D2 = DI2_D2 + DD_HS2*HS2_D2 + DD_US2*US2_D2
       DI2_MS = DI2_MS + DD_HS2*HS2_MS + DD_US2*US2_MS
       DI2_RE = DI2_RE + DD_HS2*HS2_RE + DD_US2*US2_RE
C
C----- add laminar stress contribution to outer layer CD
c###
       DD     =  0.15*(0.995-US2)**2 / RT2  * 2.0/HS2
       DD_US2 = -0.15*(0.995-US2)*2. / RT2  * 2.0/HS2
       DD_HS2 = -DD/HS2
       DD_RT2 = -DD/RT2
C
       DI2    = DI2    + DD
       DI2_U2 = DI2_U2 + DD_HS2*HS2_U2 + DD_US2*US2_U2 + DD_RT2*RT2_U2
       DI2_T2 = DI2_T2 + DD_HS2*HS2_T2 + DD_US2*US2_T2 + DD_RT2*RT2_T2
       DI2_D2 = DI2_D2 + DD_HS2*HS2_D2 + DD_US2*US2_D2
       DI2_MS = DI2_MS + DD_HS2*HS2_MS + DD_US2*US2_MS + DD_RT2*RT2_MS
       DI2_RE = DI2_RE + DD_HS2*HS2_RE + DD_US2*US2_RE + DD_RT2*RT2_RE
C
      ENDIF
C
C
      IF(ITYP.EQ.2) THEN
        CALL DIL( HK2, RT2, DI2L, DI2L_HK2, DI2L_RT2 )
C
        IF(DI2L.GT.DI2) THEN
C------- laminar CD is greater than turbulent CD -- use laminar
C-       (this will only occur for unreasonably small Rtheta)
ccc       write(*,*) 'CDt CDl Rt Hk:', DI2, DI2L, RT2, HK2
          DI2    = DI2L
          DI2_S2 = 0.
          DI2_U2 = DI2L_HK2*HK2_U2 + DI2L_RT2*RT2_U2
          DI2_T2 = DI2L_HK2*HK2_T2 + DI2L_RT2*RT2_T2
          DI2_D2 = DI2L_HK2*HK2_D2
          DI2_MS = DI2L_HK2*HK2_MS + DI2L_RT2*RT2_MS
          DI2_RE =                   DI2L_RT2*RT2_RE
        ENDIF
      ENDIF
C
cC----- add on CD contribution of inner shear layer
c       IF(ITYP.EQ.3 .AND. DW2.GT.0.0) THEN
c        DKON = 0.03*0.75**3
c        DDI = DKON*US2**3
c        DDI_US2 = 3.0*DKON*US2**2
c        DI2 = DI2 + DDI * DW2/DWTE
c        DI2_U2 = DI2_U2 + DDI_US2*US2_U2 * DW2/DWTE
c        DI2_T2 = DI2_T2 + DDI_US2*US2_T2 * DW2/DWTE
c        DI2_D2 = DI2_D2 + DDI_US2*US2_D2 * DW2/DWTE
c        DI2_MS = DI2_MS + DDI_US2*US2_MS * DW2/DWTE
c        DI2_RE = DI2_RE + DDI_US2*US2_RE * DW2/DWTE
c       ENDIF
C
      IF(ITYP.EQ.3) THEN
C------ laminar wake CD
        CALL DILW( HK2, RT2, DI2L, DI2L_HK2, DI2L_RT2 )
        IF(DI2L .GT. DI2) THEN
C------- laminar wake CD is greater than turbulent CD -- use laminar
C-       (this will only occur for unreasonably small Rtheta)
ccc         write(*,*) 'CDt CDl Rt Hk:', DI2, DI2L, RT2, HK2
         DI2    = DI2L
         DI2_S2 = 0.
         DI2_U2 = DI2L_HK2*HK2_U2 + DI2L_RT2*RT2_U2
         DI2_T2 = DI2L_HK2*HK2_T2 + DI2L_RT2*RT2_T2
         DI2_D2 = DI2L_HK2*HK2_D2
         DI2_MS = DI2L_HK2*HK2_MS + DI2L_RT2*RT2_MS
         DI2_RE =                   DI2L_RT2*RT2_RE
        ENDIF
      ENDIF
C
C
      IF(ITYP.EQ.3) THEN
C----- double dissipation for the wake (two wake halves)
       DI2    = DI2   *2.0
       DI2_S2 = DI2_S2*2.0
       DI2_U2 = DI2_U2*2.0
       DI2_T2 = DI2_T2*2.0
       DI2_D2 = DI2_D2*2.0
       DI2_MS = DI2_MS*2.0
       DI2_RE = DI2_RE*2.0
      ENDIF
C
C---- BL thickness (Delta) from simplified Green's correlation
      DE2     = (3.15 + 1.72/(HK2-1.0)   )*T2  +  D2
      DE2_HK2 = (     - 1.72/(HK2-1.0)**2)*T2
C
      DE2_U2 = DE2_HK2*HK2_U2
      DE2_T2 = DE2_HK2*HK2_T2 + (3.15 + 1.72/(HK2-1.0))
      DE2_D2 = DE2_HK2*HK2_D2 + 1.0
      DE2_MS = DE2_HK2*HK2_MS
C
ccc      HDMAX = 15.0
      HDMAX = 12.0
      IF(DE2 .GT. HDMAX*T2) THEN
cccc      IF(DE2 .GT. HDMAX*T2 .AND. (HK2 .GT. 4.0 .OR. ITYP.EQ.3)) THEN
       DE2    = HDMAX*T2
       DE2_U2 =  0.0
       DE2_T2 = HDMAX
       DE2_D2 =  0.0
       DE2_MS =  0.0
      ENDIF
C
      RETURN
      END

C===================================================================70
C
C     Sets up "dummy" BL system between airfoil TE point 
C     and first wake point infinitesimally behind TE.
C
C===================================================================70
      SUBROUTINE TESYS(CTE,TTE,DTE)

      use blpar_inc
      use xbl_inc
      IMPLICIT REAL (M)
C
      DO 55 K=1, 4
        VSREZ(K) = 0.
        VSM(K)   = 0.
        VSR(K)   = 0.
        VSX(K)   = 0.
        DO 551 L=1, 5
          VS1(K,L) = 0.
          VS2(K,L) = 0.
  551   CONTINUE
   55 CONTINUE
C
      CALL BLVAR(3)
C
      VS1(1,1) = -1.0
      VS2(1,1) = 1.0
      VSREZ(1) = CTE - S2      
C
      VS1(2,2) = -1.0
      VS2(2,2) = 1.0
      VSREZ(2) = TTE - T2
C
      VS1(3,3) = -1.0
      VS2(3,3) = 1.0
      VSREZ(3) = DTE - D2 - DW2
C
      RETURN
      END

C===================================================================70
C
C     Calculates midpoint skin friction CFM
C
C      ITYP = 1 :  laminar
C      ITYP = 2 :  turbulent
C      ITYP = 3 :  turbulent wake
C
C===================================================================70
      SUBROUTINE BLMID(ITYP)

      use blpar_inc
      use xbl_inc
      IMPLICIT REAL(M)
C
C---- set similarity variables if not defined
      IF(SIMI) THEN
       HK1    = HK2
       HK1_T1 = HK2_T2
       HK1_D1 = HK2_D2
       HK1_U1 = HK2_U2
       HK1_MS = HK2_MS
       RT1    = RT2
       RT1_T1 = RT2_T2
       RT1_U1 = RT2_U2
       RT1_MS = RT2_MS
       RT1_RE = RT2_RE
       M1    = M2
       M1_U1 = M2_U2
       M1_MS = M2_MS
      ENDIF
C
C---- define stuff for midpoint CF
      HKA = 0.5*(HK1 + HK2)
      RTA = 0.5*(RT1 + RT2)
      MA  = 0.5*(M1  + M2 )
C
C---- midpoint skin friction coefficient  (zero in wake)
      IF(ITYP.EQ.3) THEN
       CFM     = 0.
       CFM_HKA = 0.
       CFM_RTA = 0.
       CFM_MA  = 0.
       CFM_MS  = 0.
      ELSE IF(ITYP.EQ.1) THEN
       CALL CFL( HKA, RTA, MA, CFM, CFM_HKA, CFM_RTA, CFM_MA )
      ELSE
       CALL CFT( HKA, RTA, MA, CFM, CFM_HKA, CFM_RTA, CFM_MA )
       CALL CFL( HKA, RTA, MA, CFML,CFML_HKA,CFML_RTA,CFML_MA)
       IF(CFML.GT.CFM) THEN
ccc      write(*,*) 'Cft Cfl Rt Hk:', CFM, CFML, RTA, HKA, 0.5*(X1+X2)
         CFM     = CFML
         CFM_HKA = CFML_HKA
         CFM_RTA = CFML_RTA
         CFM_MA  = CFML_MA
       ENDIF
      ENDIF
C
      CFM_U1 = 0.5*(CFM_HKA*HK1_U1 + CFM_MA*M1_U1 + CFM_RTA*RT1_U1)
      CFM_T1 = 0.5*(CFM_HKA*HK1_T1 +                CFM_RTA*RT1_T1)
      CFM_D1 = 0.5*(CFM_HKA*HK1_D1                                )
C
      CFM_U2 = 0.5*(CFM_HKA*HK2_U2 + CFM_MA*M2_U2 + CFM_RTA*RT2_U2)
      CFM_T2 = 0.5*(CFM_HKA*HK2_T2 +                CFM_RTA*RT2_T2)
      CFM_D2 = 0.5*(CFM_HKA*HK2_D2                                )
C
      CFM_MS = 0.5*(CFM_HKA*HK1_MS + CFM_MA*M1_MS + CFM_RTA*RT1_MS
     &            + CFM_HKA*HK2_MS + CFM_MA*M2_MS + CFM_RTA*RT2_MS)
      CFM_RE = 0.5*(                                CFM_RTA*RT1_RE
     &                                            + CFM_RTA*RT2_RE)
C
      RETURN
      END ! BLMID

C===================================================================70
C
C     Sets up the Newton system coefficients and residuals
C
C        ITYP = 0 :  similarity station
C        ITYP = 1 :  laminar interval
C        ITYP = 2 :  turbulent interval
C        ITYP = 3 :  wake interval
C
C     This routine knows nothing about a transition interval,
C     which is taken care of by TRDIF.
C
C===================================================================70
      SUBROUTINE BLDIF(ITYP)

      use blpar_inc
      use xbl_inc
      IMPLICIT REAL(M)
C
      IF(ITYP.EQ.0) THEN
C----- similarity logarithmic differences  (prescribed)
       XLOG = 1.0
       ULOG = BULE
       TLOG = 0.5*(1.0 - BULE)
       HLOG = 0.
       DDLOG = 0.
      ELSE
C----- usual logarithmic differences
       XLOG = LOG(X2/X1)
       ULOG = LOG(U2/U1)
       TLOG = LOG(T2/T1)
       HLOG = LOG(HS2/HS1)
C       XLOG = 2.0*(X2-X1)/(X2+X1)
C       ULOG = 2.0*(U2-U1)/(U2+U1)
C       TLOG = 2.0*(T2-T1)/(T2+T1)
C       HLOG = 2.0*(HS2-HS1)/(HS2+HS1)
       DDLOG = 1.0
      ENDIF
C
      DO 55 K=1, 4
        VSREZ(K) = 0.
        VSM(K) = 0.
        VSR(K) = 0.
        VSX(K) = 0.
        DO 551 L=1, 5
          VS1(K,L) = 0.
          VS2(K,L) = 0.
  551   CONTINUE
   55 CONTINUE
C
C---- set triggering constant for local upwinding
      HUPWT = 1.0
C
ccc      HDCON = 5.0*HUPWT
ccc      HD_HK1 = 0.0
ccc      HD_HK2 = 0.0
C
      HDCON  =  5.0*HUPWT/HK2**2
      HD_HK1 =  0.0
      HD_HK2 = -HDCON*2.0/HK2
C
C---- use less upwinding in the wake
      IF(ITYP.EQ.3) THEN
       HDCON  =  HUPWT/HK2**2
       HD_HK1 =  0.0
       HD_HK2 = -HDCON*2.0/HK2
      ENDIF
C
C---- local upwinding is based on local change in  log(Hk-1)
C-    (mainly kicks in at transition)
      ARG = ABS((HK2-1.0)/(HK1-1.0))
      HL = LOG(ARG)
      HL_HK1 = -1.0/(HK1-1.0)
      HL_HK2 =  1.0/(HK2-1.0)
C
C---- set local upwinding parameter UPW and linearize it
C
C       UPW = 0.5   Trapezoidal
C       UPW = 1.0   Backward Euler
C
      HLSQ = MIN( HL**2 , 15.0 )
      EHH = EXP(-HLSQ*HDCON)
      UPW = 1.0 - 0.5*EHH
      UPW_HL =        EHH * HL  *HDCON
      UPW_HD =    0.5*EHH * HLSQ
C
      UPW_HK1 = UPW_HL*HL_HK1 + UPW_HD*HD_HK1
      UPW_HK2 = UPW_HL*HL_HK2 + UPW_HD*HD_HK2
C
      UPW_U1 = UPW_HK1*HK1_U1
      UPW_T1 = UPW_HK1*HK1_T1
      UPW_D1 = UPW_HK1*HK1_D1
      UPW_U2 = UPW_HK2*HK2_U2
      UPW_T2 = UPW_HK2*HK2_T2
      UPW_D2 = UPW_HK2*HK2_D2
      UPW_MS = UPW_HK1*HK1_MS
     &       + UPW_HK2*HK2_MS
C
C
      IF(ITYP.EQ.0) THEN
C
C***** LE point -->  set zero amplification factor
       VS2(1,1) = 1.0
       VSR(1)   = 0.
       VSREZ(1) = -AMPL2
C
      ELSE IF(ITYP.EQ.1) THEN
C
C***** laminar part -->  set amplification equation
C
C----- set average amplification AX over interval X1..X2
       CALL AXSET( HK1,    T1,    RT1, AMPL1,  
     &             HK2,    T2,    RT2, AMPL2, AMCRIT, IDAMPV,
     &      AX, AX_HK1, AX_T1, AX_RT1, AX_A1,
     &          AX_HK2, AX_T2, AX_RT2, AX_A2 )
C
       REZC = AMPL2 - AMPL1 - AX*(X2-X1)
       Z_AX = -(X2-X1)
C
       VS1(1,1) = Z_AX* AX_A1  -  1.0
       VS1(1,2) = Z_AX*(AX_HK1*HK1_T1 + AX_T1 + AX_RT1*RT1_T1)
       VS1(1,3) = Z_AX*(AX_HK1*HK1_D1                        )
       VS1(1,4) = Z_AX*(AX_HK1*HK1_U1         + AX_RT1*RT1_U1)
       VS1(1,5) =  AX
       VS2(1,1) = Z_AX* AX_A2  +  1.0
       VS2(1,2) = Z_AX*(AX_HK2*HK2_T2 + AX_T2 + AX_RT2*RT2_T2)
       VS2(1,3) = Z_AX*(AX_HK2*HK2_D2                        )         
       VS2(1,4) = Z_AX*(AX_HK2*HK2_U2         + AX_RT2*RT2_U2)
       VS2(1,5) = -AX
       VSM(1)   = Z_AX*(AX_HK1*HK1_MS         + AX_RT1*RT1_MS
     &                + AX_HK2*HK2_MS         + AX_RT2*RT2_MS)
       VSR(1)   = Z_AX*(                        AX_RT1*RT1_RE
     &                                        + AX_RT2*RT2_RE)
       VSX(1)   = 0.
       VSREZ(1) = -REZC
C
      ELSE
C
C***** turbulent part -->  set shear lag equation
C
       SA  = (1.0-UPW)*S1  + UPW*S2
       CQA = (1.0-UPW)*CQ1 + UPW*CQ2
       CFA = (1.0-UPW)*CF1 + UPW*CF2
       HKA = (1.0-UPW)*HK1 + UPW*HK2
C
       USA = 0.5*(US1 + US2)
       RTA = 0.5*(RT1 + RT2)
       DEA = 0.5*(DE1 + DE2)
       DA  = 0.5*(D1  + D2 )
C
C
       IF(ITYP.EQ.3) THEN
C------ increased dissipation length in wake (decrease its reciprocal)
        ALD = DLCON
       ELSE
        ALD = 1.0
       ENDIF
C
C----- set and linearize  equilibrium 1/Ue dUe/dx   ...  NEW  12 Oct 94
       IF(ITYP.EQ.2) THEN
        GCC = GCCON
        HKC     = HKA - 1.0 - GCC/RTA
        HKC_HKA = 1.0
        HKC_RTA =             GCC/RTA**2
        IF(HKC .LT. 0.01) THEN
         HKC = 0.01
         HKC_HKA = 0.0
         HKC_RTA = 0.0
        ENDIF
       ELSE
        GCC = 0.0
        HKC = HKA - 1.0
        HKC_HKA = 1.0
        HKC_RTA = 0.0
       ENDIF
C
       HR     = HKC     / (GACON*ALD*HKA)
       HR_HKA = HKC_HKA / (GACON*ALD*HKA) - HR / HKA
       HR_RTA = HKC_RTA / (GACON*ALD*HKA)
C
       UQ     = (0.5*CFA - HR**2) / (GBCON*DA)
       UQ_HKA =   -2.0*HR*HR_HKA  / (GBCON*DA)
       UQ_RTA =   -2.0*HR*HR_RTA  / (GBCON*DA)
       UQ_CFA =   0.5             / (GBCON*DA)
       UQ_DA  = -UQ/DA
       UQ_UPW = UQ_CFA*(CF2-CF1) + UQ_HKA*(HK2-HK1)
C
       UQ_T1 = (1.0-UPW)*(UQ_CFA*CF1_T1 + UQ_HKA*HK1_T1) + UQ_UPW*UPW_T1
       UQ_D1 = (1.0-UPW)*(UQ_CFA*CF1_D1 + UQ_HKA*HK1_D1) + UQ_UPW*UPW_D1
       UQ_U1 = (1.0-UPW)*(UQ_CFA*CF1_U1 + UQ_HKA*HK1_U1) + UQ_UPW*UPW_U1
       UQ_T2 =      UPW *(UQ_CFA*CF2_T2 + UQ_HKA*HK2_T2) + UQ_UPW*UPW_T2
       UQ_D2 =      UPW *(UQ_CFA*CF2_D2 + UQ_HKA*HK2_D2) + UQ_UPW*UPW_D2
       UQ_U2 =      UPW *(UQ_CFA*CF2_U2 + UQ_HKA*HK2_U2) + UQ_UPW*UPW_U2
       UQ_MS = (1.0-UPW)*(UQ_CFA*CF1_MS + UQ_HKA*HK1_MS) + UQ_UPW*UPW_MS
     &       +      UPW *(UQ_CFA*CF2_MS + UQ_HKA*HK2_MS)
       UQ_RE = (1.0-UPW)* UQ_CFA*CF1_RE
     &       +      UPW * UQ_CFA*CF2_RE
C
       UQ_T1 = UQ_T1             + 0.5*UQ_RTA*RT1_T1
       UQ_D1 = UQ_D1 + 0.5*UQ_DA
       UQ_U1 = UQ_U1             + 0.5*UQ_RTA*RT1_U1
       UQ_T2 = UQ_T2             + 0.5*UQ_RTA*RT2_T2
       UQ_D2 = UQ_D2 + 0.5*UQ_DA
       UQ_U2 = UQ_U2             + 0.5*UQ_RTA*RT2_U2
       UQ_MS = UQ_MS             + 0.5*UQ_RTA*RT1_MS
     &                           + 0.5*UQ_RTA*RT2_MS
       UQ_RE = UQ_RE             + 0.5*UQ_RTA*RT1_RE
     &                           + 0.5*UQ_RTA*RT2_RE
C
       SCC = SCCON*1.333/(1.0+USA)
       SCC_USA = -SCC/(1.0+USA)
C
       SCC_US1 = SCC_USA*0.5
       SCC_US2 = SCC_USA*0.5
C
C
       SLOG = LOG(S2/S1)
       DXI = X2 - X1
C
       REZC = SCC*(CQA - SA*ALD)*DXI
     &      - DEA*2.0*          SLOG
     &      + DEA*2.0*(UQ*DXI - ULOG)*DUXCON
C

c        if(  ! (rt2.gt.1.0e3 .and. rt1.le.1.0e3) .or.
c     &     (rt2.gt.1.0e4 .and. rt1.le.1.0e4) .or.
c     &     (rt2.gt.1.0e5 .and. rt1.le.1.0e5)        ) then
c           gga = (HKA-1.0-GCC/RTA)/HKA / sqrt(0.5*CFA)
c           write(*,4455) rta, hka, gga, cfa, cqa, sa, uq, ulog/dxi
c 4455      format(1x,f7.0, 2f9.4,f10.6,2f8.5,2f10.5)
c        endif


       Z_CFA = DEA*2.0*UQ_CFA*DXI * DUXCON
       Z_HKA = DEA*2.0*UQ_HKA*DXI * DUXCON
       Z_DA  = DEA*2.0*UQ_DA *DXI * DUXCON
       Z_SL = -DEA*2.0
       Z_UL = -DEA*2.0 * DUXCON
       Z_DXI = SCC    *(CQA - SA*ALD)     + DEA*2.0*UQ*DUXCON
       Z_USA = SCC_USA*(CQA - SA*ALD)*DXI
       Z_CQA = SCC*DXI
       Z_SA = -SCC*DXI*ALD
       Z_DEA = 2.0*((UQ*DXI - ULOG)*DUXCON - SLOG)
C
       Z_UPW = Z_CQA*(CQ2-CQ1) + Z_SA *(S2 -S1 )
     &       + Z_CFA*(CF2-CF1) + Z_HKA*(HK2-HK1)
       Z_DE1 = 0.5*Z_DEA
       Z_DE2 = 0.5*Z_DEA
       Z_US1 = 0.5*Z_USA
       Z_US2 = 0.5*Z_USA
       Z_D1  = 0.5*Z_DA
       Z_D2  = 0.5*Z_DA
       Z_U1  =                 - Z_UL/U1
       Z_U2  =                   Z_UL/U2
       Z_X1  = -Z_DXI
       Z_X2  =  Z_DXI
       Z_S1  = (1.0-UPW)*Z_SA  - Z_SL/S1
       Z_S2  =      UPW *Z_SA  + Z_SL/S2
       Z_CQ1 = (1.0-UPW)*Z_CQA
       Z_CQ2 =      UPW *Z_CQA
       Z_CF1 = (1.0-UPW)*Z_CFA
       Z_CF2 =      UPW *Z_CFA
       Z_HK1 = (1.0-UPW)*Z_HKA
       Z_HK2 =      UPW *Z_HKA
C
       VS1(1,1) = Z_S1
       VS1(1,2) =        Z_UPW*UPW_T1 + Z_DE1*DE1_T1 + Z_US1*US1_T1
       VS1(1,3) = Z_D1 + Z_UPW*UPW_D1 + Z_DE1*DE1_D1 + Z_US1*US1_D1
       VS1(1,4) = Z_U1 + Z_UPW*UPW_U1 + Z_DE1*DE1_U1 + Z_US1*US1_U1
       VS1(1,5) = Z_X1
       VS2(1,1) = Z_S2
       VS2(1,2) =        Z_UPW*UPW_T2 + Z_DE2*DE2_T2 + Z_US2*US2_T2
       VS2(1,3) = Z_D2 + Z_UPW*UPW_D2 + Z_DE2*DE2_D2 + Z_US2*US2_D2
       VS2(1,4) = Z_U2 + Z_UPW*UPW_U2 + Z_DE2*DE2_U2 + Z_US2*US2_U2
       VS2(1,5) = Z_X2
       VSM(1)   =        Z_UPW*UPW_MS + Z_DE1*DE1_MS + Z_US1*US1_MS
     &                                + Z_DE2*DE2_MS + Z_US2*US2_MS
C
       VS1(1,2) = VS1(1,2) + Z_CQ1*CQ1_T1 + Z_CF1*CF1_T1 + Z_HK1*HK1_T1
       VS1(1,3) = VS1(1,3) + Z_CQ1*CQ1_D1 + Z_CF1*CF1_D1 + Z_HK1*HK1_D1
       VS1(1,4) = VS1(1,4) + Z_CQ1*CQ1_U1 + Z_CF1*CF1_U1 + Z_HK1*HK1_U1
C
       VS2(1,2) = VS2(1,2) + Z_CQ2*CQ2_T2 + Z_CF2*CF2_T2 + Z_HK2*HK2_T2
       VS2(1,3) = VS2(1,3) + Z_CQ2*CQ2_D2 + Z_CF2*CF2_D2 + Z_HK2*HK2_D2
       VS2(1,4) = VS2(1,4) + Z_CQ2*CQ2_U2 + Z_CF2*CF2_U2 + Z_HK2*HK2_U2
C
       VSM(1)   = VSM(1)   + Z_CQ1*CQ1_MS + Z_CF1*CF1_MS + Z_HK1*HK1_MS
     &                     + Z_CQ2*CQ2_MS + Z_CF2*CF2_MS + Z_HK2*HK2_MS
       VSR(1)   =            Z_CQ1*CQ1_RE + Z_CF1*CF1_RE
     &                     + Z_CQ2*CQ2_RE + Z_CF2*CF2_RE
       VSX(1)   = 0.
       VSREZ(1) = -REZC
C
      ENDIF
C
C**** Set up momentum equation
      HA = 0.5*(H1 + H2)
      MA = 0.5*(M1 + M2)
      XA = 0.5*(X1 + X2)
      TA = 0.5*(T1 + T2)
      HWA = 0.5*(DW1/T1 + DW2/T2)
C
C---- set Cf term, using central value CFM for better accuracy in drag
      CFX     = 0.50*CFM*XA/TA  +  0.25*(CF1*X1/T1 + CF2*X2/T2)
      CFX_XA  = 0.50*CFM   /TA
      CFX_TA  = -.50*CFM*XA/TA**2
C
      CFX_X1  = 0.25*CF1   /T1     + CFX_XA*0.5
      CFX_X2  = 0.25*CF2   /T2     + CFX_XA*0.5
      CFX_T1  = -.25*CF1*X1/T1**2  + CFX_TA*0.5
      CFX_T2  = -.25*CF2*X2/T2**2  + CFX_TA*0.5
      CFX_CF1 = 0.25*    X1/T1
      CFX_CF2 = 0.25*    X2/T2
      CFX_CFM = 0.50*    XA/TA
C
      BTMP = HA + 2.0 - MA + HWA
C
      REZT  = TLOG + BTMP*ULOG - XLOG*0.5*CFX
      Z_CFX = -XLOG*0.5
      Z_HA  =  ULOG
      Z_HWA =  ULOG
      Z_MA  = -ULOG
      Z_XL  =-DDLOG * 0.5*CFX
      Z_UL  = DDLOG * BTMP
      Z_TL  = DDLOG
C
      Z_CFM = Z_CFX*CFX_CFM
      Z_CF1 = Z_CFX*CFX_CF1
      Z_CF2 = Z_CFX*CFX_CF2
C
      Z_T1 = -Z_TL/T1 + Z_CFX*CFX_T1 + Z_HWA*0.5*(-DW1/T1**2)
      Z_T2 =  Z_TL/T2 + Z_CFX*CFX_T2 + Z_HWA*0.5*(-DW2/T2**2)
      Z_X1 = -Z_XL/X1 + Z_CFX*CFX_X1
      Z_X2 =  Z_XL/X2 + Z_CFX*CFX_X2
      Z_U1 = -Z_UL/U1
      Z_U2 =  Z_UL/U2
C
      VS1(2,2) = 0.5*Z_HA*H1_T1 + Z_CFM*CFM_T1 + Z_CF1*CF1_T1 + Z_T1
      VS1(2,3) = 0.5*Z_HA*H1_D1 + Z_CFM*CFM_D1 + Z_CF1*CF1_D1
      VS1(2,4) = 0.5*Z_MA*M1_U1 + Z_CFM*CFM_U1 + Z_CF1*CF1_U1 + Z_U1
      VS1(2,5) =                                                Z_X1
      VS2(2,2) = 0.5*Z_HA*H2_T2 + Z_CFM*CFM_T2 + Z_CF2*CF2_T2 + Z_T2
      VS2(2,3) = 0.5*Z_HA*H2_D2 + Z_CFM*CFM_D2 + Z_CF2*CF2_D2
      VS2(2,4) = 0.5*Z_MA*M2_U2 + Z_CFM*CFM_U2 + Z_CF2*CF2_U2 + Z_U2
      VS2(2,5) =                                                Z_X2
C
      VSM(2)   = 0.5*Z_MA*M1_MS + Z_CFM*CFM_MS + Z_CF1*CF1_MS
     &         + 0.5*Z_MA*M2_MS                + Z_CF2*CF2_MS
      VSR(2)   =                  Z_CFM*CFM_RE + Z_CF1*CF1_RE
     &                                         + Z_CF2*CF2_RE
      VSX(2)   = 0.
      VSREZ(2) = -REZT
C
C**** Set up shape parameter equation
C
      XOT1 = X1/T1
      XOT2 = X2/T2
C
      HA  = 0.5*(H1  + H2 )
      HSA = 0.5*(HS1 + HS2)
      HCA = 0.5*(HC1 + HC2)
      HWA = 0.5*(DW1/T1 + DW2/T2)
C
      DIX = (1.0-UPW)*DI1*XOT1 + UPW*DI2*XOT2
      CFX = (1.0-UPW)*CF1*XOT1 + UPW*CF2*XOT2
      DIX_UPW = DI2*XOT2 - DI1*XOT1
      CFX_UPW = CF2*XOT2 - CF1*XOT1
C
      BTMP = 2.0*HCA/HSA + 1.0 - HA - HWA
C
      REZH  = HLOG + BTMP*ULOG + XLOG*(0.5*CFX-DIX)
      Z_CFX =  XLOG*0.5
      Z_DIX = -XLOG
      Z_HCA = 2.0*ULOG/HSA
      Z_HA  = -ULOG
      Z_HWA = -ULOG
      Z_XL  = DDLOG * (0.5*CFX-DIX)
      Z_UL  = DDLOG * BTMP
      Z_HL  = DDLOG
C
      Z_UPW = Z_CFX*CFX_UPW + Z_DIX*DIX_UPW
C
      Z_HS1 = -HCA*ULOG/HSA**2 - Z_HL/HS1
      Z_HS2 = -HCA*ULOG/HSA**2 + Z_HL/HS2
C
      Z_CF1 = (1.0-UPW)*Z_CFX*XOT1
      Z_CF2 =      UPW *Z_CFX*XOT2
      Z_DI1 = (1.0-UPW)*Z_DIX*XOT1
      Z_DI2 =      UPW *Z_DIX*XOT2
C
      Z_T1 = (1.0-UPW)*(Z_CFX*CF1 + Z_DIX*DI1)*(-XOT1/T1)
      Z_T2 =      UPW *(Z_CFX*CF2 + Z_DIX*DI2)*(-XOT2/T2)
      Z_X1 = (1.0-UPW)*(Z_CFX*CF1 + Z_DIX*DI1)/ T1        - Z_XL/X1
      Z_X2 =      UPW *(Z_CFX*CF2 + Z_DIX*DI2)/ T2        + Z_XL/X2
      Z_U1 =                                              - Z_UL/U1
      Z_U2 =                                                Z_UL/U2
C
      Z_T1 = Z_T1 + Z_HWA*0.5*(-DW1/T1**2)
      Z_T2 = Z_T2 + Z_HWA*0.5*(-DW2/T2**2)
C
      VS1(3,1) =                               Z_DI1*DI1_S1
      VS1(3,2) = Z_HS1*HS1_T1 + Z_CF1*CF1_T1 + Z_DI1*DI1_T1 + Z_T1
      VS1(3,3) = Z_HS1*HS1_D1 + Z_CF1*CF1_D1 + Z_DI1*DI1_D1
      VS1(3,4) = Z_HS1*HS1_U1 + Z_CF1*CF1_U1 + Z_DI1*DI1_U1 + Z_U1
      VS1(3,5) =                                              Z_X1
      VS2(3,1) =                               Z_DI2*DI2_S2
      VS2(3,2) = Z_HS2*HS2_T2 + Z_CF2*CF2_T2 + Z_DI2*DI2_T2 + Z_T2
      VS2(3,3) = Z_HS2*HS2_D2 + Z_CF2*CF2_D2 + Z_DI2*DI2_D2
      VS2(3,4) = Z_HS2*HS2_U2 + Z_CF2*CF2_U2 + Z_DI2*DI2_U2 + Z_U2
      VS2(3,5) =                                              Z_X2
      VSM(3)   = Z_HS1*HS1_MS + Z_CF1*CF1_MS + Z_DI1*DI1_MS
     &         + Z_HS2*HS2_MS + Z_CF2*CF2_MS + Z_DI2*DI2_MS
      VSR(3)   = Z_HS1*HS1_RE + Z_CF1*CF1_RE + Z_DI1*DI1_RE
     &         + Z_HS2*HS2_RE + Z_CF2*CF2_RE + Z_DI2*DI2_RE
C
      VS1(3,2) = VS1(3,2) + 0.5*(Z_HCA*HC1_T1+Z_HA*H1_T1) + Z_UPW*UPW_T1
      VS1(3,3) = VS1(3,3) + 0.5*(Z_HCA*HC1_D1+Z_HA*H1_D1) + Z_UPW*UPW_D1
      VS1(3,4) = VS1(3,4) + 0.5*(Z_HCA*HC1_U1           ) + Z_UPW*UPW_U1
      VS2(3,2) = VS2(3,2) + 0.5*(Z_HCA*HC2_T2+Z_HA*H2_T2) + Z_UPW*UPW_T2
      VS2(3,3) = VS2(3,3) + 0.5*(Z_HCA*HC2_D2+Z_HA*H2_D2) + Z_UPW*UPW_D2
      VS2(3,4) = VS2(3,4) + 0.5*(Z_HCA*HC2_U2           ) + Z_UPW*UPW_U2
C
      VSM(3)   = VSM(3)   + 0.5*(Z_HCA*HC1_MS           ) + Z_UPW*UPW_MS
     &                    + 0.5*(Z_HCA*HC2_MS           )
C
      VSX(3)   = 0.
      VSREZ(3) = -REZH
C
      RETURN
      END

C===================================================================70
C
C     Sets up the Newton system governing the
C     transition interval.  Equations governing
C     the  laminar  part  X1 < xi < XT  and
C     the turbulent part  XT < xi < X2
C     are simply summed.
C
C===================================================================70
      SUBROUTINE TRDIF

      use blpar_inc
      use xbl_inc
      IMPLICIT REAL(M)
      REAL*8 :: BL1(4,5), BL2(4,5), BLREZ(4), BLM(4), BLR(4), BLX(4)
     &    , BT1(4,5), BT2(4,5), BTREZ(4), BTM(4), BTR(4), BTX(4)
C
C---- save variables and sensitivities for future restoration
C     DP mod: to remove need for EQUIVALENCE and COM1, COM2
      call store_c1sav()
      call store_c2sav()
C      DO 5 ICOM=1, NCOM
C        C1SAV(ICOM) = COM1(ICOM)
C        C2SAV(ICOM) = COM2(ICOM)
C    5 CONTINUE
C
C---- weighting factors for linear interpolation to transition point
      WF2    = (XT-X1)/(X2-X1)
      WF2_XT = 1.0/(X2-X1)
C
      WF2_A1 = WF2_XT*XT_A1
      WF2_X1 = WF2_XT*XT_X1 + (WF2-1.0)/(X2-X1)
      WF2_X2 = WF2_XT*XT_X2 -  WF2     /(X2-X1)
      WF2_T1 = WF2_XT*XT_T1
      WF2_T2 = WF2_XT*XT_T2
      WF2_D1 = WF2_XT*XT_D1
      WF2_D2 = WF2_XT*XT_D2
      WF2_U1 = WF2_XT*XT_U1
      WF2_U2 = WF2_XT*XT_U2
      WF2_MS = WF2_XT*XT_MS
      WF2_RE = WF2_XT*XT_RE
      WF2_XF = WF2_XT*XT_XF
C
      WF1    = 1.0 - WF2
      WF1_A1 = -WF2_A1
      WF1_X1 = -WF2_X1
      WF1_X2 = -WF2_X2
      WF1_T1 = -WF2_T1
      WF1_T2 = -WF2_T2
      WF1_D1 = -WF2_D1
      WF1_D2 = -WF2_D2
      WF1_U1 = -WF2_U1
      WF1_U2 = -WF2_U2
      WF1_MS = -WF2_MS
      WF1_RE = -WF2_RE
      WF1_XF = -WF2_XF
C
C
C**** FIRST,  do laminar part between X1 and XT
C
C-----interpolate primary variables to transition point
      TT    = T1*WF1    + T2*WF2
      TT_A1 = T1*WF1_A1 + T2*WF2_A1
      TT_X1 = T1*WF1_X1 + T2*WF2_X1
      TT_X2 = T1*WF1_X2 + T2*WF2_X2
      TT_T1 = T1*WF1_T1 + T2*WF2_T1 + WF1
      TT_T2 = T1*WF1_T2 + T2*WF2_T2 + WF2
      TT_D1 = T1*WF1_D1 + T2*WF2_D1
      TT_D2 = T1*WF1_D2 + T2*WF2_D2
      TT_U1 = T1*WF1_U1 + T2*WF2_U1
      TT_U2 = T1*WF1_U2 + T2*WF2_U2
      TT_MS = T1*WF1_MS + T2*WF2_MS
      TT_RE = T1*WF1_RE + T2*WF2_RE
      TT_XF = T1*WF1_XF + T2*WF2_XF
C
      DT    = D1*WF1    + D2*WF2
      DT_A1 = D1*WF1_A1 + D2*WF2_A1
      DT_X1 = D1*WF1_X1 + D2*WF2_X1
      DT_X2 = D1*WF1_X2 + D2*WF2_X2
      DT_T1 = D1*WF1_T1 + D2*WF2_T1
      DT_T2 = D1*WF1_T2 + D2*WF2_T2
      DT_D1 = D1*WF1_D1 + D2*WF2_D1 + WF1
      DT_D2 = D1*WF1_D2 + D2*WF2_D2 + WF2
      DT_U1 = D1*WF1_U1 + D2*WF2_U1
      DT_U2 = D1*WF1_U2 + D2*WF2_U2
      DT_MS = D1*WF1_MS + D2*WF2_MS
      DT_RE = D1*WF1_RE + D2*WF2_RE
      DT_XF = D1*WF1_XF + D2*WF2_XF
C
      UT    = U1*WF1    + U2*WF2
      UT_A1 = U1*WF1_A1 + U2*WF2_A1
      UT_X1 = U1*WF1_X1 + U2*WF2_X1
      UT_X2 = U1*WF1_X2 + U2*WF2_X2
      UT_T1 = U1*WF1_T1 + U2*WF2_T1
      UT_T2 = U1*WF1_T2 + U2*WF2_T2
      UT_D1 = U1*WF1_D1 + U2*WF2_D1
      UT_D2 = U1*WF1_D2 + U2*WF2_D2
      UT_U1 = U1*WF1_U1 + U2*WF2_U1 + WF1
      UT_U2 = U1*WF1_U2 + U2*WF2_U2 + WF2
      UT_MS = U1*WF1_MS + U2*WF2_MS
      UT_RE = U1*WF1_RE + U2*WF2_RE
      UT_XF = U1*WF1_XF + U2*WF2_XF
C
C---- set primary "T" variables at XT  (really placed into "2" variables)
      X2 = XT
      T2 = TT
      D2 = DT
      U2 = UT
C
      AMPL2 = AMCRIT
      S2 = 0.
C
C---- calculate laminar secondary "T" variables
      CALL BLKIN
      CALL BLVAR(1)
C
C---- calculate X1-XT midpoint CFM value
      CALL BLMID(1)
C=
C=    at this point, all "2" variables are really "T" variables at XT
C=
C
C---- set up Newton system for dAm, dTh, dDs, dUe, dXi  at  X1 and XT
      CALL BLDIF(1)
C
C---- The current Newton system is in terms of "1" and "T" variables,
C-    so calculate its equivalent in terms of "1" and "2" variables.
C-    In other words, convert residual sensitivities wrt "T" variables
C-    into sensitivities wrt "1" and "2" variables.  The amplification
C-    equation is unnecessary here, so the K=1 row is left empty.
      DO 10 K=2, 3
        BLREZ(K) = VSREZ(K)
        BLM(K)   = VSM(K)
     &           + VS2(K,2)*TT_MS
     &           + VS2(K,3)*DT_MS
     &           + VS2(K,4)*UT_MS
     &           + VS2(K,5)*XT_MS
        BLR(K)   = VSR(K)
     &           + VS2(K,2)*TT_RE
     &           + VS2(K,3)*DT_RE
     &           + VS2(K,4)*UT_RE
     &           + VS2(K,5)*XT_RE
        BLX(K)   = VSX(K)
     &           + VS2(K,2)*TT_XF
     &           + VS2(K,3)*DT_XF
     &           + VS2(K,4)*UT_XF
     &           + VS2(K,5)*XT_XF
C
        BL1(K,1) = VS1(K,1)
     &           + VS2(K,2)*TT_A1
     &           + VS2(K,3)*DT_A1
     &           + VS2(K,4)*UT_A1
     &           + VS2(K,5)*XT_A1
        BL1(K,2) = VS1(K,2)
     &           + VS2(K,2)*TT_T1
     &           + VS2(K,3)*DT_T1
     &           + VS2(K,4)*UT_T1
     &           + VS2(K,5)*XT_T1
        BL1(K,3) = VS1(K,3)
     &           + VS2(K,2)*TT_D1
     &           + VS2(K,3)*DT_D1
     &           + VS2(K,4)*UT_D1
     &           + VS2(K,5)*XT_D1
        BL1(K,4) = VS1(K,4)
     &           + VS2(K,2)*TT_U1
     &           + VS2(K,3)*DT_U1
     &           + VS2(K,4)*UT_U1
     &           + VS2(K,5)*XT_U1
        BL1(K,5) = VS1(K,5)
     &           + VS2(K,2)*TT_X1
     &           + VS2(K,3)*DT_X1
     &           + VS2(K,4)*UT_X1
     &           + VS2(K,5)*XT_X1
C
        BL2(K,1) = 0.
        BL2(K,2) = VS2(K,2)*TT_T2
     &           + VS2(K,3)*DT_T2
     &           + VS2(K,4)*UT_T2
     &           + VS2(K,5)*XT_T2
        BL2(K,3) = VS2(K,2)*TT_D2
     &           + VS2(K,3)*DT_D2
     &           + VS2(K,4)*UT_D2
     &           + VS2(K,5)*XT_D2
        BL2(K,4) = VS2(K,2)*TT_U2
     &           + VS2(K,3)*DT_U2
     &           + VS2(K,4)*UT_U2
     &           + VS2(K,5)*XT_U2
        BL2(K,5) = VS2(K,2)*TT_X2
     &           + VS2(K,3)*DT_X2
     &           + VS2(K,4)*UT_X2
     &           + VS2(K,5)*XT_X2
C
   10 CONTINUE
C
C
C**** SECOND, set up turbulent part between XT and X2  ****
C
C---- calculate equilibrium shear coefficient CQT at transition point
      CALL BLVAR(2)
C
C---- set initial shear coefficient value ST at transition point
C-    ( note that CQ2, CQ2_T2, etc. are really "CQT", "CQT_TT", etc.)
C
      CTR     = CTRCON*EXP(-CTRCEX/(HK2-1.0))
      CTR_HK2 = CTR * CTRCEX/(HK2-1.0)**2
C
c      CTR     = 1.1*EXP(-10.0/HK2**2)
c      CTR_HK2 = CTR * 10.0 * 2.0/HK2**3
C
CCC      CTR = 1.2
CCC      CTR = 0.7
CCC      CTR_HK2 = 0.0
C
      ST    = CTR*CQ2
      ST_TT = CTR*CQ2_T2 + CQ2*CTR_HK2*HK2_T2
      ST_DT = CTR*CQ2_D2 + CQ2*CTR_HK2*HK2_D2
      ST_UT = CTR*CQ2_U2 + CQ2*CTR_HK2*HK2_U2
      ST_MS = CTR*CQ2_MS + CQ2*CTR_HK2*HK2_MS
      ST_RE = CTR*CQ2_RE
C
C---- calculate ST sensitivities wrt the actual "1" and "2" variables
      ST_A1 = ST_TT*TT_A1 + ST_DT*DT_A1 + ST_UT*UT_A1
      ST_X1 = ST_TT*TT_X1 + ST_DT*DT_X1 + ST_UT*UT_X1
      ST_X2 = ST_TT*TT_X2 + ST_DT*DT_X2 + ST_UT*UT_X2
      ST_T1 = ST_TT*TT_T1 + ST_DT*DT_T1 + ST_UT*UT_T1
      ST_T2 = ST_TT*TT_T2 + ST_DT*DT_T2 + ST_UT*UT_T2
      ST_D1 = ST_TT*TT_D1 + ST_DT*DT_D1 + ST_UT*UT_D1
      ST_D2 = ST_TT*TT_D2 + ST_DT*DT_D2 + ST_UT*UT_D2
      ST_U1 = ST_TT*TT_U1 + ST_DT*DT_U1 + ST_UT*UT_U1
      ST_U2 = ST_TT*TT_U2 + ST_DT*DT_U2 + ST_UT*UT_U2
      ST_MS = ST_TT*TT_MS + ST_DT*DT_MS + ST_UT*UT_MS + ST_MS
      ST_RE = ST_TT*TT_RE + ST_DT*DT_RE + ST_UT*UT_RE + ST_RE
      ST_XF = ST_TT*TT_XF + ST_DT*DT_XF + ST_UT*UT_XF
C
      AMPL2 = 0.
      S2 = ST
C
C---- recalculate turbulent secondary "T" variables using proper CTI
      CALL BLVAR(2)
C
C---- set "1" variables to "T" variables and reset "2" variables
C-    to their saved turbulent values
C
C     DP mod: to remove need for EQUIVALENCE COM1, COM2, and COMMON
      call com2_to_com1()
      call from_c2sav()
C      DO 30 ICOM=1, NCOM
C        COM1(ICOM) = COM2(ICOM)
C        COM2(ICOM) = C2SAV(ICOM)
C   30 CONTINUE
C
C---- calculate XT-X2 midpoint CFM value
      CALL BLMID(2)
C
C---- set up Newton system for dCt, dTh, dDs, dUe, dXi  at  XT and X2
      CALL BLDIF(2)
C
C---- convert sensitivities wrt "T" variables into sensitivities
C-    wrt "1" and "2" variables as done before for the laminar part
      DO 40 K=1, 3
        BTREZ(K) = VSREZ(K)
        BTM(K)   = VSM(K) 
     &           + VS1(K,1)*ST_MS
     &           + VS1(K,2)*TT_MS
     &           + VS1(K,3)*DT_MS
     &           + VS1(K,4)*UT_MS
     &           + VS1(K,5)*XT_MS
        BTR(K)   = VSR(K) 
     &           + VS1(K,1)*ST_RE
     &           + VS1(K,2)*TT_RE
     &           + VS1(K,3)*DT_RE
     &           + VS1(K,4)*UT_RE
     &           + VS1(K,5)*XT_RE
        BTX(K)   = VSX(K)
     &           + VS1(K,1)*ST_XF
     &           + VS1(K,2)*TT_XF
     &           + VS1(K,3)*DT_XF
     &           + VS1(K,4)*UT_XF
     &           + VS1(K,5)*XT_XF
C
        BT1(K,1) = VS1(K,1)*ST_A1
     &           + VS1(K,2)*TT_A1
     &           + VS1(K,3)*DT_A1
     &           + VS1(K,4)*UT_A1
     &           + VS1(K,5)*XT_A1
        BT1(K,2) = VS1(K,1)*ST_T1
     &           + VS1(K,2)*TT_T1
     &           + VS1(K,3)*DT_T1
     &           + VS1(K,4)*UT_T1
     &           + VS1(K,5)*XT_T1
        BT1(K,3) = VS1(K,1)*ST_D1
     &           + VS1(K,2)*TT_D1
     &           + VS1(K,3)*DT_D1
     &           + VS1(K,4)*UT_D1
     &           + VS1(K,5)*XT_D1
        BT1(K,4) = VS1(K,1)*ST_U1
     &           + VS1(K,2)*TT_U1
     &           + VS1(K,3)*DT_U1
     &           + VS1(K,4)*UT_U1
     &           + VS1(K,5)*XT_U1
        BT1(K,5) = VS1(K,1)*ST_X1
     &           + VS1(K,2)*TT_X1
     &           + VS1(K,3)*DT_X1
     &           + VS1(K,4)*UT_X1
     &           + VS1(K,5)*XT_X1
C
        BT2(K,1) = VS2(K,1)
        BT2(K,2) = VS2(K,2)
     &           + VS1(K,1)*ST_T2
     &           + VS1(K,2)*TT_T2
     &           + VS1(K,3)*DT_T2
     &           + VS1(K,4)*UT_T2
     &           + VS1(K,5)*XT_T2
        BT2(K,3) = VS2(K,3)
     &           + VS1(K,1)*ST_D2
     &           + VS1(K,2)*TT_D2
     &           + VS1(K,3)*DT_D2
     &           + VS1(K,4)*UT_D2
     &           + VS1(K,5)*XT_D2
        BT2(K,4) = VS2(K,4)
     &           + VS1(K,1)*ST_U2
     &           + VS1(K,2)*TT_U2
     &           + VS1(K,3)*DT_U2
     &           + VS1(K,4)*UT_U2
     &           + VS1(K,5)*XT_U2
        BT2(K,5) = VS2(K,5)
     &           + VS1(K,1)*ST_X2
     &           + VS1(K,2)*TT_X2
     &           + VS1(K,3)*DT_X2
     &           + VS1(K,4)*UT_X2
     &           + VS1(K,5)*XT_X2
C
   40 CONTINUE
C
C---- Add up laminar and turbulent parts to get final system
C-    in terms of honest-to-God "1" and "2" variables.
      VSREZ(1) =            BTREZ(1)
      VSREZ(2) = BLREZ(2) + BTREZ(2)
      VSREZ(3) = BLREZ(3) + BTREZ(3)
      VSM(1)   =            BTM(1)
      VSM(2)   = BLM(2)   + BTM(2)
      VSM(3)   = BLM(3)   + BTM(3)
      VSR(1)   =            BTR(1)
      VSR(2)   = BLR(2)   + BTR(2)
      VSR(3)   = BLR(3)   + BTR(3)
      VSX(1)   =            BTX(1)
      VSX(2)   = BLX(2)   + BTX(2)
      VSX(3)   = BLX(3)   + BTX(3)
      DO 60 L=1, 5
        VS1(1,L) =            BT1(1,L)
        VS2(1,L) =            BT2(1,L)
        VS1(2,L) = BL1(2,L) + BT1(2,L)
        VS2(2,L) = BL2(2,L) + BT2(2,L)
        VS1(3,L) = BL1(3,L) + BT1(3,L)
        VS2(3,L) = BL2(3,L) + BT2(3,L)
   60 CONTINUE
C
C---- To be sanitary, restore "1" quantities which got clobbered
C-    in all of the numerical gymnastics above.  The "2" variables
C-    were already restored for the XT-X2 differencing part.
C     DP mod: to remove need for EQUIVALENCE and COM1
      call from_c1sav()
C      DO 70 ICOM=1, NCOM
C        COM1(ICOM) = C1SAV(ICOM)
C   70 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C
C     Sets up the BL Newton system governing the current interval:
C
C     |       ||dA1|     |       ||dA2|       |     |
C     |  VS1  ||dT1|  +  |  VS2  ||dT2|   =   |VSREZ|
C     |       ||dD1|     |       ||dD2|       |     |
C              |dU1|              |dU2|
C              |dX1|              |dX2|
C
C        3x5    5x1         3x5    5x1          3x1
C
C     The system as shown corresponds to a laminar station
C     If TRAN, then  dS2  replaces  dA2
C     If TURB, then  dS1, dS2  replace  dA1, dA2
C
C
C===================================================================70
      SUBROUTINE BLSYS

      use blpar_inc
      use xbl_inc
      IMPLICIT REAL(M)
C
C---- calculate secondary BL variables and their sensitivities
      IF(WAKE) THEN
       CALL BLVAR(3)
       CALL BLMID(3)
      ELSE IF(TURB.OR.TRAN) THEN
       CALL BLVAR(2)
       CALL BLMID(2)
      ELSE
       CALL BLVAR(1)
       CALL BLMID(1)
      ENDIF
C
C---- for the similarity station, "1" and "2" variables are the same
      IF(SIMI) THEN
C      DP mod: to remove need for EQUIVALENCE COM1, COM2, and COMMON
       call com2_to_com1()
C       DO 3 ICOM=1, NCOM
C         COM1(ICOM) = COM2(ICOM)
C    3  CONTINUE
      ENDIF
C
C---- set up appropriate finite difference system for current interval
      IF(TRAN) THEN
       CALL TRDIF
      ELSE IF(SIMI) THEN
       CALL BLDIF(0)
      ELSE IF(.NOT.TURB) THEN
       CALL BLDIF(1)
      ELSE IF(WAKE) THEN
       CALL BLDIF(3)
      ELSE IF(TURB) THEN
       CALL BLDIF(2)
      ENDIF
C
      IF(SIMI) THEN
C----- at similarity station, "1" variables are really "2" variables
       DO 10 K=1, 4
         DO 101 L=1, 5
           VS2(K,L) = VS1(K,L) + VS2(K,L)
           VS1(K,L) = 0.
  101    CONTINUE
   10  CONTINUE
      ENDIF
C
C---- change system over into incompressible Uei and Mach
      DO 20 K=1, 4
C
C------ residual derivatives wrt compressible Uec
        RES_U1 = VS1(K,4)
        RES_U2 = VS2(K,4)
        RES_MS = VSM(K)
C
C------ combine with derivatives of compressible  U1,U2 = Uec(Uei M)
        VS1(K,4) = RES_U1*U1_UEI
        VS2(K,4) =                RES_U2*U2_UEI
        VSM(K)   = RES_U1*U1_MS + RES_U2*U2_MS  + RES_MS
   20 CONTINUE
C
      RETURN
      END

C===================================================================70
C===================================================================70
      SUBROUTINE DSLIM(DSTR,THET,UEDG,MSQ,HKLIM)
      IMPLICIT REAL (A-H,M,O-Z)
C
      H = DSTR/THET
      CALL HKIN(H,MSQ,HK,HK_H,HK_M)
C
      DH = MAX( 0.0 , HKLIM-HK ) / HK_H
      DSTR = DSTR + DH*THET
C
      RETURN
      END

C===================================================================70
C
C     Marches the BLs and wake in direct mode using
C     the UEDG array. If separation is encountered,
C     a plausible value of Hk extrapolated from
C     upstream is prescribed instead.  Continuous
C     checking of transition onset is performed.
C
C===================================================================70
      SUBROUTINE MRCHUE

      use xfoil_inc
      use blpar_inc
      use xbl_inc

      LOGICAL :: DIRECT
      REAL*8 :: MSQ
C
C---- shape parameters for separation criteria
      HLMAX = 3.8
      HTMAX = 2.5
C
      DO 2000 IS=1, 2
C
C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE) WRITE(*,*) '   side ', IS, ' ...'
C
C---- set forced transition arc length position
      CALL XIFSET(IS)
C
C---- initialize similarity station with Thwaites' formula
      IBL = 2
      XSI = XSSI(IBL,IS)
      UEI = UEDG(IBL,IS)
C      BULE = LOG(UEDG(IBL+1,IS)/UEI) / LOG(XSSI(IBL+1,IS)/XSI)
C      BULE = MAX( -.08 , BULE )
      BULE = 1.0
      UCON = UEI/XSI**BULE
      TSQ = 0.45/(UCON*(5.0*BULE+1.0)*REYBL) * XSI**(1.0-BULE)
      THI = SQRT(TSQ)
      DSI = 2.2*THI
      AMI = 0.0
C
C---- initialize Ctau for first turbulent station
      CTI = 0.03
C
      TRAN = .FALSE.
      TURB = .FALSE.
      ITRAN(IS) = IBLTE(IS)
C
C---- march downstream
      DO 1000 IBL=2, NBL(IS)
        IBM = IBL-1
C
        IW = IBL - IBLTE(IS)
C
        SIMI = IBL.EQ.2
        WAKE = IBL.GT.IBLTE(IS)
C
C------ prescribed quantities
        XSI = XSSI(IBL,IS)
        UEI = UEDG(IBL,IS)
C
        IF(WAKE) THEN
         IW = IBL - IBLTE(IS)
         DSWAKI = WGAP(IW)
        ELSE
         DSWAKI = 0.
        ENDIF
C
        DIRECT = .TRUE.
C
C------ Newton iteration loop for current station
        DO 100 ITBL=1, 25
C
C-------- assemble 10x3 linearized system for dCtau, dTh, dDs, dUe, dXi
C         at the previous "1" station and the current "2" station
C         (the "1" station coefficients will be ignored)
C
C
          CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
          CALL BLKIN
C
C-------- check for transition and set appropriate flags and things
          IF((.NOT.SIMI) .AND. (.NOT.TURB)) THEN
           CALL TRCHEK
C          DP mod: check for infinite loop condition
           IF (XFOIL_FAIL) RETURN
           AMI = AMPL2
C
           IF(TRAN) THEN
            ITRAN(IS) = IBL
            IF(CTI.LE.0.0) THEN
             CTI = 0.03
             S2 = CTI
            ENDIF
           ELSE
            ITRAN(IS) = IBL+2
           ENDIF
C
C
          ENDIF
C
          IF(IBL.EQ.IBLTE(IS)+1) THEN
           TTE = THET(IBLTE(1),1) + THET(IBLTE(2),2)
           DTE = DSTR(IBLTE(1),1) + DSTR(IBLTE(2),2) + ANTE
           CTE = ( CTAU(IBLTE(1),1)*THET(IBLTE(1),1)
     &           + CTAU(IBLTE(2),2)*THET(IBLTE(2),2) ) / TTE
           CALL TESYS(CTE,TTE,DTE)
          ELSE
           CALL BLSYS
          ENDIF
C
          IF(DIRECT) THEN
C
C--------- try direct mode (set dUe = 0 in currently empty 4th line)
           VS2(4,1) = 0.
           VS2(4,2) = 0.
           VS2(4,3) = 0.
           VS2(4,4) = 1.0
           VSREZ(4) = 0.
C
C--------- solve Newton system for current "2" station
           CALL GAUSS(4,4,VS2,VSREZ,1)
C
C--------- determine max changes and underrelax if necessary
           DMAX = MAX( ABS(VSREZ(2)/THI),
     &                 ABS(VSREZ(3)/DSI)  )
           IF(IBL.LT.ITRAN(IS)) DMAX = MAX(DMAX,ABS(VSREZ(1)/10.0))
           IF(IBL.GE.ITRAN(IS)) DMAX = MAX(DMAX,ABS(VSREZ(1)/CTI ))
C
           RLX = 1.0
           IF(DMAX.GT.0.3) RLX = 0.3/DMAX
C
C--------- see if direct mode is not applicable
           IF(IBL .NE. IBLTE(IS)+1) THEN
C
C---------- calculate resulting kinematic shape parameter Hk
            MSQ = UEI*UEI*HSTINV / (GM1BL*(1.0 - 0.5*UEI*UEI*HSTINV))
            HTEST = (DSI + RLX*VSREZ(3)) / (THI + RLX*VSREZ(2))
            CALL HKIN( HTEST, MSQ, HKTEST, DUMMY, DUMMY)
C
C---------- decide whether to do direct or inverse problem based on Hk
            IF(IBL.LT.ITRAN(IS)) HMAX = HLMAX
            IF(IBL.GE.ITRAN(IS)) HMAX = HTMAX
            DIRECT = HKTEST.LT.HMAX
           ENDIF
C
           IF(DIRECT) THEN
C---------- update as usual
ccc            IF(IBL.LT.ITRAN(IS)) AMI = AMI + RLX*VSREZ(1)
            IF(IBL.GE.ITRAN(IS)) CTI = CTI + RLX*VSREZ(1)
            THI = THI + RLX*VSREZ(2)
            DSI = DSI + RLX*VSREZ(3)
           ELSE
C---------- set prescribed Hk for inverse calculation at the current station
            IF(IBL.LT.ITRAN(IS)) THEN
C----------- laminar case: relatively slow increase in Hk downstream
             HTARG = HK1 + 0.03*(X2-X1)/T1
            ELSE IF(IBL.EQ.ITRAN(IS)) THEN
C----------- transition interval: weighted laminar and turbulent case
             HTARG = HK1 + (0.03*(XT-X1) - 0.15*(X2-XT))/T1
            ELSE IF(WAKE) THEN
C----------- turbulent wake case:
C-           asymptotic wake behavior with approximate Backward Euler
             CONST = 0.03*(X2-X1)/T1
             HK2 = HK1
             HK2 = HK2 - (HK2 +     CONST*(HK2-1.0)**3 - HK1)
     &                  /(1.0 + 3.0*CONST*(HK2-1.0)**2)
             HK2 = HK2 - (HK2 +     CONST*(HK2-1.0)**3 - HK1)
     &                  /(1.0 + 3.0*CONST*(HK2-1.0)**2)
             HK2 = HK2 - (HK2 +     CONST*(HK2-1.0)**3 - HK1)
     &                  /(1.0 + 3.0*CONST*(HK2-1.0)**2)
             HTARG = HK2
            ELSE
C----------- turbulent case: relatively fast decrease in Hk downstream
             HTARG = HK1 - 0.15*(X2-X1)/T1
            ENDIF
C
C---------- limit specified Hk to something reasonable
            IF(WAKE) THEN
             HTARG = MAX( HTARG , 1.01 )
            ELSE
             HTARG = MAX( HTARG , HMAX )
            ENDIF
C
C           DP mod: added SILENT_MODE option
            IF (.NOT. SILENT_MODE) WRITE(*,1300) IBL, HTARG
 1300       FORMAT(' MRCHUE: Inverse mode at', I4, '     Hk =', F8.3)
C
C---------- try again with prescribed Hk
            GO TO 100
C
           ENDIF
C
          ELSE
C
C-------- inverse mode (force Hk to prescribed value HTARG)
           VS2(4,1) = 0.
           VS2(4,2) = HK2_T2
           VS2(4,3) = HK2_D2
           VS2(4,4) = HK2_U2
           VSREZ(4) = HTARG - HK2
C
           CALL GAUSS(4,4,VS2,VSREZ,1)
C
C--------- added Ue clamp   MD  3 Apr 03
           DMAX = MAX( ABS(VSREZ(2)/THI),
     &                 ABS(VSREZ(3)/DSI),
     &                 ABS(VSREZ(4)/UEI)  )
           IF(IBL.GE.ITRAN(IS)) DMAX = MAX( DMAX , ABS(VSREZ(1)/CTI))
C
           RLX = 1.0
           IF(DMAX.GT.0.3) RLX = 0.3/DMAX
C
C--------- update variables
ccc           IF(IBL.LT.ITRAN(IS)) AMI = AMI + RLX*VSREZ(1)
           IF(IBL.GE.ITRAN(IS)) CTI = CTI + RLX*VSREZ(1)
           THI = THI + RLX*VSREZ(2)
           DSI = DSI + RLX*VSREZ(3)
           UEI = UEI + RLX*VSREZ(4)
C
          ENDIF
C
C-------- eliminate absurd transients
          IF(IBL.GE.ITRAN(IS)) THEN
           CTI = MIN(CTI , 0.30 )
           CTI = MAX(CTI , 0.0000001 )
          ENDIF
C
          IF(IBL.LE.IBLTE(IS)) THEN
            HKLIM = 1.02
          ELSE
            HKLIM = 1.00005
          ENDIF
          MSQ = UEI*UEI*HSTINV / (GM1BL*(1.0 - 0.5*UEI*UEI*HSTINV))
          DSW = DSI - DSWAKI
          CALL DSLIM(DSW,THI,UEI,MSQ,HKLIM)
          DSI = DSW + DSWAKI
C
          IF(DMAX.LE.1.0E-5) GO TO 110
C
  100   CONTINUE
C       DP mod: added SILENT_MODE option
        IF (.NOT. SILENT_MODE) WRITE(*,1350) IBL, IS, DMAX 
 1350   FORMAT(' MRCHUE: Convergence failed at',I4,'  side',I2,
     &         '    Res =', E12.4)
C
C------ the current unconverged solution might still be reasonable...
CCC        IF(DMAX .LE. 0.1) GO TO 110
        IF(DMAX .LE. 0.1) GO TO 109
C
C------- the current solution is garbage --> extrapolate values instead
         IF(IBL.GT.3) THEN 
          IF(IBL.LE.IBLTE(IS)) THEN
           THI = THET(IBM,IS) * (XSSI(IBL,IS)/XSSI(IBM,IS))**0.5
           DSI = DSTR(IBM,IS) * (XSSI(IBL,IS)/XSSI(IBM,IS))**0.5
          ELSE IF(IBL.EQ.IBLTE(IS)+1) THEN
           CTI = CTE
           THI = TTE
           DSI = DTE
          ELSE
           THI = THET(IBM,IS)
           RATLEN = (XSSI(IBL,IS)-XSSI(IBM,IS)) / (10.0*DSTR(IBM,IS))
           DSI = (DSTR(IBM,IS) + THI*RATLEN) / (1.0 + RATLEN)
          ENDIF
          IF(IBL.EQ.ITRAN(IS)) CTI = 0.05
          IF(IBL.GT.ITRAN(IS)) CTI = CTAU(IBM,IS)
C
          UEI = UEDG(IBL,IS)
          IF(IBL.GT.2 .AND. IBL.LT.NBL(IS))
     &     UEI = 0.5*(UEDG(IBL-1,IS) + UEDG(IBL+1,IS))
         ENDIF
C
 109     CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
         CALL BLKIN
C
C------- check for transition and set appropriate flags and things
         IF((.NOT.SIMI) .AND. (.NOT.TURB)) THEN
          CALL TRCHEK
C         DP mod: check for infinite loop condition
          IF (XFOIL_FAIL) RETURN
          AMI = AMPL2
          IF(     TRAN) ITRAN(IS) = IBL
          IF(.NOT.TRAN) ITRAN(IS) = IBL+2
         ENDIF
C
C------- set all other extrapolated values for current station
         IF(IBL.LT.ITRAN(IS)) CALL BLVAR(1)
         IF(IBL.GE.ITRAN(IS)) CALL BLVAR(2)
         IF(WAKE) CALL BLVAR(3)
C
         IF(IBL.LT.ITRAN(IS)) CALL BLMID(1)
         IF(IBL.GE.ITRAN(IS)) CALL BLMID(2)
         IF(WAKE) CALL BLMID(3)
C
C------ pick up here after the Newton iterations
  110   CONTINUE
C
C------ store primary variables
        IF(IBL.LT.ITRAN(IS)) CTAU(IBL,IS) = AMI
        IF(IBL.GE.ITRAN(IS)) CTAU(IBL,IS) = CTI
        THET(IBL,IS) = THI
        DSTR(IBL,IS) = DSI
        UEDG(IBL,IS) = UEI
        MASS(IBL,IS) = DSI*UEI
        TAU(IBL,IS)  = 0.5*R2*U2*U2*CF2
        DIS(IBL,IS)  =     R2*U2*U2*U2*DI2*HS2*0.5
        CTQ(IBL,IS)  = CQ2
        DELT(IBL,IS) = DE2
        TSTR(IBL,IS) = HS2*T2
C
C------ set "1" variables to "2" variables for next streamwise station
        CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
        CALL BLKIN
C       DP mod: to remove need for EQUIVALENCE COM1, COM2, and COMMON
        call com2_to_com1()
C        DO 310 ICOM=1, NCOM
C          COM1(ICOM) = COM2(ICOM)
C  310   CONTINUE
C
C------ turbulent intervals will follow transition interval or TE
        IF(TRAN .OR. IBL.EQ.IBLTE(IS)) THEN
         TURB = .TRUE.
C
C------- save transition location
         TFORCE(IS) = TRFORC
         XSSITR(IS) = XT
        ENDIF
C
        TRAN = .FALSE.
C
        IF(IBL.EQ.IBLTE(IS)) THEN
         THI = THET(IBLTE(1),1) + THET(IBLTE(2),2)
         DSI = DSTR(IBLTE(1),1) + DSTR(IBLTE(2),2) + ANTE
        ENDIF
C
 1000 CONTINUE
 2000 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Marches the BLs and wake in mixed mode using
C     the current Ue and Hk.  The calculated Ue
C     and Hk lie along a line quasi-normal to the
C     natural Ue-Hk characteristic line of the
C     current BL so that the Goldstein or Levy-Lees
C     singularity is never encountered.  Continuous
C     checking of transition onset is performed.
C
C===================================================================70
      SUBROUTINE MRCHDU

      use xfoil_inc
      use blpar_inc
      use xbl_inc

      REAL*8 :: VTMP(4,5), VZTMP(4)
      REAL*8 :: MSQ
ccc   REAL MDI
C
      DATA DEPS / 5.0E-6 /
C
C---- constant controlling how far Hk is allowed to deviate
C-    from the specified value.
      SENSWT = 1000.0
C
      DO 2000 IS=1, 2
C
C---- set forced transition arc length position
      CALL XIFSET(IS)
C
C---- set leading edge pressure gradient parameter  x/u du/dx
      IBL = 2
      XSI = XSSI(IBL,IS)
      UEI = UEDG(IBL,IS)
CCC      BULE = LOG(UEDG(IBL+1,IS)/UEI) / LOG(XSSI(IBL+1,IS)/XSI)
CCC      BULE = MAX( -.08 , BULE )
      BULE = 1.0
C
C---- old transition station
      ITROLD = ITRAN(IS)
C
      TRAN = .FALSE.
      TURB = .FALSE.
      ITRAN(IS) = IBLTE(IS)
C
C---- march downstream
      DO 1000 IBL=2, NBL(IS)
        IBM = IBL-1
C
        SIMI = IBL.EQ.2
        WAKE = IBL.GT.IBLTE(IS)
C
C------ initialize current station to existing variables
        XSI = XSSI(IBL,IS)
        UEI = UEDG(IBL,IS)
        THI = THET(IBL,IS)
        DSI = DSTR(IBL,IS)

CCC        MDI = MASS(IBL,IS)
C
C------ fixed BUG   MD 7 June 99
        IF(IBL.LT.ITROLD) THEN
         AMI = CTAU(IBL,IS)
         CTI = 0.03
        ELSE
         CTI = CTAU(IBL,IS)
         IF(CTI.LE.0.0) CTI = 0.03
        ENDIF
C
CCC        DSI = MDI/UEI
C
        IF(WAKE) THEN
         IW = IBL - IBLTE(IS)
         DSWAKI = WGAP(IW)
        ELSE
         DSWAKI = 0.
        ENDIF
C
        IF(IBL.LE.IBLTE(IS)) DSI = MAX(DSI-DSWAKI,1.02000*THI) + DSWAKI
        IF(IBL.GT.IBLTE(IS)) DSI = MAX(DSI-DSWAKI,1.00005*THI) + DSWAKI
C
C------ Newton iteration loop for current station
        DO 100 ITBL=1, 25
C
C-------- assemble 10x3 linearized system for dCtau, dTh, dDs, dUe, dXi
C         at the previous "1" station and the current "2" station
C         (the "1" station coefficients will be ignored)
C
          CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
          CALL BLKIN
C
C-------- check for transition and set appropriate flags and things
          IF((.NOT.SIMI) .AND. (.NOT.TURB)) THEN
           CALL TRCHEK
C          DP mod: check for infinite loop condition
           IF (XFOIL_FAIL) RETURN
           AMI = AMPL2
           IF(     TRAN) ITRAN(IS) = IBL
           IF(.NOT.TRAN) ITRAN(IS) = IBL+2
          ENDIF
C
          IF(IBL.EQ.IBLTE(IS)+1) THEN
           TTE = THET(IBLTE(1),1) + THET(IBLTE(2),2)
           DTE = DSTR(IBLTE(1),1) + DSTR(IBLTE(2),2) + ANTE
           CTE = ( CTAU(IBLTE(1),1)*THET(IBLTE(1),1)
     &           + CTAU(IBLTE(2),2)*THET(IBLTE(2),2) ) / TTE
           CALL TESYS(CTE,TTE,DTE)
          ELSE
           CALL BLSYS
          ENDIF
C
C-------- set stuff at first iteration...
          IF(ITBL.EQ.1) THEN
C
C--------- set "baseline" Ue and Hk for forming  Ue(Hk)  relation
           UEREF = U2
           HKREF = HK2
C
C--------- if current point IBL was turbulent and is now laminar, then...
           IF(IBL.LT.ITRAN(IS) .AND. IBL.GE.ITROLD ) THEN
C---------- extrapolate baseline Hk
            UEM = UEDG(IBL-1,IS)
            DSM = DSTR(IBL-1,IS)
            THM = THET(IBL-1,IS)
            MSQ = UEM*UEM*HSTINV / (GM1BL*(1.0 - 0.5*UEM*UEM*HSTINV))
            CALL HKIN( DSM/THM, MSQ, HKREF, DUMMY, DUMMY )
           ENDIF
C
C--------- if current point IBL was laminar, then...
           IF(IBL.LT.ITROLD) THEN
C---------- reinitialize or extrapolate Ctau if it's now turbulent
            IF(TRAN) CTAU(IBL,IS) = 0.03
            IF(TURB) CTAU(IBL,IS) = CTAU(IBL-1,IS)
            IF(TRAN .OR. TURB) THEN
             CTI = CTAU(IBL,IS)
             S2 = CTI
            ENDIF
           ENDIF
C
          ENDIF
C
C
          IF(SIMI .OR. IBL.EQ.IBLTE(IS)+1) THEN
C
C--------- for similarity station or first wake point, prescribe Ue
           VS2(4,1) = 0.
           VS2(4,2) = 0.
           VS2(4,3) = 0.
           VS2(4,4) = U2_UEI
           VSREZ(4) = UEREF - U2
C
          ELSE
C
C********* calculate Ue-Hk characteristic slope
C
           DO 20 K=1, 4
             VZTMP(K) = VSREZ(K)
             DO 201 L=1, 5
               VTMP(K,L) = VS2(K,L)
  201        CONTINUE
   20      CONTINUE
C
C--------- set unit dHk
           VTMP(4,1) = 0.
           VTMP(4,2) = HK2_T2
           VTMP(4,3) = HK2_D2
           VTMP(4,4) = HK2_U2*U2_UEI
           VZTMP(4)  = 1.0
C
C--------- calculate dUe response
           CALL GAUSS(4,4,VTMP,VZTMP,1)
C
C--------- set  SENSWT * (normalized dUe/dHk)
           SENNEW = SENSWT * VZTMP(4) * HKREF/UEREF
           IF(ITBL.LE.5) THEN
            SENS = SENNEW
           ELSE IF(ITBL.LE.15) THEN
            SENS = 0.5*(SENS + SENNEW)
           ENDIF
C
C--------- set prescribed Ue-Hk combination
           VS2(4,1) = 0.
           VS2(4,2) =  HK2_T2 * HKREF
           VS2(4,3) =  HK2_D2 * HKREF
           VS2(4,4) =( HK2_U2 * HKREF  +  SENS/UEREF )*U2_UEI
           VSREZ(4) = -(HKREF**2)*(HK2 / HKREF - 1.0)
     &                     - SENS*(U2  / UEREF - 1.0)
C
          ENDIF
C
C-------- solve Newton system for current "2" station
          CALL GAUSS(4,4,VS2,VSREZ,1)
C
C-------- determine max changes and underrelax if necessary
C-------- (added Ue clamp   MD  3 Apr 03)
          DMAX = MAX( ABS(VSREZ(2)/THI),
     &                ABS(VSREZ(3)/DSI),
     &                ABS(VSREZ(4)/UEI)  )
          IF(IBL.GE.ITRAN(IS)) DMAX = MAX(DMAX,ABS(VSREZ(1)/(10.0*CTI)))
C
          RLX = 1.0
          IF(DMAX.GT.0.3) RLX = 0.3/DMAX
C
C-------- update as usual
          IF(IBL.LT.ITRAN(IS)) AMI = AMI + RLX*VSREZ(1)
          IF(IBL.GE.ITRAN(IS)) CTI = CTI + RLX*VSREZ(1)
          THI = THI + RLX*VSREZ(2)
          DSI = DSI + RLX*VSREZ(3)
          UEI = UEI + RLX*VSREZ(4)
C
C-------- eliminate absurd transients
          IF(IBL.GE.ITRAN(IS)) THEN
           CTI = MIN(CTI , 0.30 )
           CTI = MAX(CTI , 0.0000001 )
          ENDIF
C
          IF(IBL.LE.IBLTE(IS)) THEN
            HKLIM = 1.02
          ELSE
            HKLIM = 1.00005
          ENDIF
          MSQ = UEI*UEI*HSTINV / (GM1BL*(1.0 - 0.5*UEI*UEI*HSTINV))
          DSW = DSI - DSWAKI
          CALL DSLIM(DSW,THI,UEI,MSQ,HKLIM)
          DSI = DSW + DSWAKI
C
          IF(DMAX.LE.DEPS) GO TO 110
C
  100   CONTINUE
C
        IF (.NOT. SILENT_MODE) WRITE(*,1350) IBL, IS, DMAX 
 1350   FORMAT(' MRCHDU: Convergence failed at',I4,'  side',I2,
     &         '    Res =', E12.4)
C
C------ the current unconverged solution might still be reasonable...
CCC        IF(DMAX .LE. 0.1) GO TO 110
        IF(DMAX .LE. 0.1) GO TO 109
C
C------- the current solution is garbage --> extrapolate values instead
         IF(IBL.GT.3) THEN
          IF(IBL.LE.IBLTE(IS)) THEN
           THI = THET(IBM,IS) * (XSSI(IBL,IS)/XSSI(IBM,IS))**0.5
           DSI = DSTR(IBM,IS) * (XSSI(IBL,IS)/XSSI(IBM,IS))**0.5
           UEI = UEDG(IBM,IS)
          ELSE IF(IBL.EQ.IBLTE(IS)+1) THEN
           CTI = CTE
           THI = TTE
           DSI = DTE
           UEI = UEDG(IBM,IS)
          ELSE
           THI = THET(IBM,IS)
           RATLEN = (XSSI(IBL,IS)-XSSI(IBM,IS)) / (10.0*DSTR(IBM,IS))
           DSI = (DSTR(IBM,IS) + THI*RATLEN) / (1.0 + RATLEN)
           UEI = UEDG(IBM,IS)
          ENDIF
          IF(IBL.EQ.ITRAN(IS)) CTI = 0.05
          IF(IBL.GT.ITRAN(IS)) CTI = CTAU(IBM,IS)
         ENDIF
C
 109     CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
         CALL BLKIN
C
C------- check for transition and set appropriate flags and things
         IF((.NOT.SIMI) .AND. (.NOT.TURB)) THEN
          CALL TRCHEK
C         DP mod: check for infinite loop condition
          IF (XFOIL_FAIL) RETURN
          AMI = AMPL2
          IF(     TRAN) ITRAN(IS) = IBL
          IF(.NOT.TRAN) ITRAN(IS) = IBL+2
         ENDIF
C
C------- set all other extrapolated values for current station
         IF(IBL.LT.ITRAN(IS)) CALL BLVAR(1)
         IF(IBL.GE.ITRAN(IS)) CALL BLVAR(2)
         IF(WAKE) CALL BLVAR(3)
C
         IF(IBL.LT.ITRAN(IS)) CALL BLMID(1)
         IF(IBL.GE.ITRAN(IS)) CALL BLMID(2)
         IF(WAKE) CALL BLMID(3)
C
C------ pick up here after the Newton iterations
  110   CONTINUE
C
        SENS = SENNEW
C
C------ store primary variables
        IF(IBL.LT.ITRAN(IS)) CTAU(IBL,IS) = AMI
        IF(IBL.GE.ITRAN(IS)) CTAU(IBL,IS) = CTI
        THET(IBL,IS) = THI
        DSTR(IBL,IS) = DSI
        UEDG(IBL,IS) = UEI
        MASS(IBL,IS) = DSI*UEI
        TAU(IBL,IS)  = 0.5*R2*U2*U2*CF2
        DIS(IBL,IS)  =     R2*U2*U2*U2*DI2*HS2*0.5
        CTQ(IBL,IS)  = CQ2
        DELT(IBL,IS) = DE2
        TSTR(IBL,IS) = HS2*T2
C
C------ set "1" variables to "2" variables for next streamwise station
        CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
        CALL BLKIN
C       DP mod: to remove need for EQUIVALENCE COM1, COM2, and COMMON
        call com2_to_com1()
C        DO 310 ICOM=1, NCOM
C          COM1(ICOM) = COM2(ICOM)
C  310   CONTINUE
C
C
C------ turbulent intervals will follow transition interval or TE
        IF(TRAN .OR. IBL.EQ.IBLTE(IS)) THEN
         TURB = .TRUE.
C
C------- save transition location
         TFORCE(IS) = TRFORC
         XSSITR(IS) = XT
        ENDIF
C
        TRAN = .FALSE.
C
 1000 CONTINUE
C
 2000 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Sets Ue from inviscid Ue plus all source influence
C
C===================================================================70
      SUBROUTINE UESET

      use xfoil_inc
C
      DO 1 IS=1, 2
        DO 10 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
C
          DUI = 0.
          DO 100 JS=1, 2
            DO 1000 JBL=2, NBL(JS)
              J  = IPAN(JBL,JS)
              UE_M = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
              DUI = DUI + UE_M*MASS(JBL,JS)
 1000       CONTINUE
  100     CONTINUE
C
          UEDG(IBL,IS) = UINV(IBL,IS) + DUI
C
   10   CONTINUE
    1 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C     Sets up the BL Newton system coefficients
C     for the current BL variables and the edge
C     velocities received from SETUP. The local
C     BL system coefficients are then
C     incorporated into the global Newton system.  
C
C===================================================================70
      SUBROUTINE SETBL

      use xfoil_inc
      use blpar_inc
      use xbl_inc

      REAL*8 :: USAV(IVX,2)
      REAL*8 :: U1_M(2*IVX), U2_M(2*IVX)
      REAL*8 :: D1_M(2*IVX), D2_M(2*IVX)
      REAL*8 :: ULE1_M(2*IVX), ULE2_M(2*IVX)
      REAL*8 :: UTE1_M(2*IVX), UTE2_M(2*IVX)
      REAL*8 :: MA_CLMR, MSQ_CLMR, MDI
      HVRAT = 0.0
C
C---- set the CL used to define Mach, Reynolds numbers
      IF(LALFA) THEN
       CLMR = CL
      ELSE
       CLMR = CLSPEC
      ENDIF
C
C---- set current MINF(CL)
      CALL MRCL(CLMR,MA_CLMR,RE_CLMR)
      MSQ_CLMR = 2.0*MINF*MA_CLMR
C
C---- set compressibility parameter TKLAM and derivative TK_MSQ
      CALL COMSET
C
C---- set gas constant (= Cp/Cv)
      GAMBL = GAMMA
      GM1BL = GAMM1
C
C---- set parameters for compressibility correction
      QINFBL = QINF
      TKBL    = TKLAM
      TKBL_MS = TKL_MSQ
C
C---- stagnation density and 1/enthalpy
      RSTBL    = (1.0 + 0.5*GM1BL*MINF**2) ** (1.0/GM1BL)
      RSTBL_MS = 0.5*RSTBL/(1.0 + 0.5*GM1BL*MINF**2)
C
      HSTINV    = GM1BL*(MINF/QINFBL)**2 / (1.0 + 0.5*GM1BL*MINF**2)
      HSTINV_MS = GM1BL*( 1.0/QINFBL)**2 / (1.0 + 0.5*GM1BL*MINF**2)
     &                - 0.5*GM1BL*HSTINV / (1.0 + 0.5*GM1BL*MINF**2)
C
C---- set Reynolds number based on freestream density, velocity, viscosity
      HERAT    = 1.0 - 0.5*QINFBL**2*HSTINV
      HERAT_MS =     - 0.5*QINFBL**2*HSTINV_MS
C
      REYBL    = REINF * SQRT(HERAT**3) * (1.0+HVRAT)/(HERAT+HVRAT)
      REYBL_RE =         SQRT(HERAT**3) * (1.0+HVRAT)/(HERAT+HVRAT)
      REYBL_MS = REYBL * (1.5/HERAT - 1.0/(HERAT+HVRAT))*HERAT_MS
C
      AMCRIT = ACRIT
!FIXME: IDAMP always set to 0 here.  In xfoil, can be set to 1 by user.
      IDAMPV = IDAMP
C
C---- save TE thickness
      DWTE = WGAP(1)
C
      IF(.NOT.LBLINI) THEN
C----- initialize BL by marching with Ue (fudge at separation)
C      DP mod: added SILENT_MODE option
       IF (.NOT. SILENT_MODE) THEN
         WRITE(*,*)
         WRITE(*,*) 'Initializing BL ...'
       ENDIF
       CALL MRCHUE
C      DP mod: check for infinite loop condition
       IF (XFOIL_FAIL) RETURN
       LBLINI = .TRUE.
      ENDIF
C
      IF (.NOT. SILENT_MODE) WRITE(*,*)
C
C---- march BL with current Ue and Ds to establish transition
      CALL MRCHDU
C     DP mod: check for infinite loop condition
      IF (XFOIL_FAIL) RETURN
C
      DO 5 IS=1, 2
        DO 6 IBL=2, NBL(IS)
          USAV(IBL,IS) = UEDG(IBL,IS)
    6   CONTINUE
    5 CONTINUE
C
      CALL UESET
C
      DO 7 IS=1, 2
        DO 8 IBL=2, NBL(IS)
          TEMP = USAV(IBL,IS)
          USAV(IBL,IS) = UEDG(IBL,IS)
          UEDG(IBL,IS) = TEMP
    8   CONTINUE
    7 CONTINUE
C
      ILE1 = IPAN(2,1)
      ILE2 = IPAN(2,2)
      ITE1 = IPAN(IBLTE(1),1)
      ITE2 = IPAN(IBLTE(2),2)
C
      JVTE1 = ISYS(IBLTE(1),1)
      JVTE2 = ISYS(IBLTE(2),2)
C
      DULE1 = UEDG(2,1) - USAV(2,1)
      DULE2 = UEDG(2,2) - USAV(2,2)
C
C---- set LE and TE Ue sensitivities wrt all m values
      DO 10 JS=1, 2
        DO 110 JBL=2, NBL(JS)
          J  = IPAN(JBL,JS)
          JV = ISYS(JBL,JS)
          ULE1_M(JV) = -VTI(       2,1)*VTI(JBL,JS)*DIJ(ILE1,J)
          ULE2_M(JV) = -VTI(       2,2)*VTI(JBL,JS)*DIJ(ILE2,J)
          UTE1_M(JV) = -VTI(IBLTE(1),1)*VTI(JBL,JS)*DIJ(ITE1,J)
          UTE2_M(JV) = -VTI(IBLTE(2),2)*VTI(JBL,JS)*DIJ(ITE2,J)
  110   CONTINUE
   10 CONTINUE
C
      ULE1_A = UINV_A(2,1)
      ULE2_A = UINV_A(2,2)
C
C**** Go over each boundary layer/wake
      DO 2000 IS=1, 2
C
C---- there is no station "1" at similarity, so zero everything out
      DO 20 JS=1, 2
        DO 210 JBL=2, NBL(JS)
          JV = ISYS(JBL,JS)
          U1_M(JV) = 0.
          D1_M(JV) = 0.
  210   CONTINUE
   20 CONTINUE
      U1_A = 0.
      D1_A = 0.
C
      DUE1 = 0.
      DDS1 = 0.
C
C---- similarity station pressure gradient parameter  x/u du/dx
      IBL = 2
      BULE = 1.0
C
C---- set forced transition arc length position
      CALL XIFSET(IS)
C
      TRAN = .FALSE.
      TURB = .FALSE.
C
C**** Sweep downstream setting up BL equation linearizations
      DO 1000 IBL=2, NBL(IS)
C
      IV  = ISYS(IBL,IS)
C
      SIMI = IBL.EQ.2
      WAKE = IBL.GT.IBLTE(IS)
      TRAN = IBL.EQ.ITRAN(IS)
      TURB = IBL.GT.ITRAN(IS)
C
      I = IPAN(IBL,IS)
C
C---- set primary variables for current station
      XSI = XSSI(IBL,IS)
      IF(IBL.LT.ITRAN(IS)) AMI = CTAU(IBL,IS)
      IF(IBL.GE.ITRAN(IS)) CTI = CTAU(IBL,IS)
      UEI = UEDG(IBL,IS)
      THI = THET(IBL,IS)
      MDI = MASS(IBL,IS)
C
      DSI = MDI/UEI
C
      IF(WAKE) THEN
       IW = IBL - IBLTE(IS)
       DSWAKI = WGAP(IW)
      ELSE
       DSWAKI = 0.
      ENDIF
C
C---- set derivatives of DSI (= D2)
      D2_M2 =  1.0/UEI
      D2_U2 = -DSI/UEI
C
      DO 30 JS=1, 2
        DO 310 JBL=2, NBL(JS)
          J  = IPAN(JBL,JS)
          JV = ISYS(JBL,JS)
          U2_M(JV) = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
          D2_M(JV) = D2_U2*U2_M(JV)
  310   CONTINUE
   30 CONTINUE
      D2_M(IV) = D2_M(IV) + D2_M2
C
      U2_A = UINV_A(IBL,IS)
      D2_A = D2_U2*U2_A
C
C---- "forced" changes due to mismatch between UEDG and USAV=UINV+dij*MASS
      DUE2 = UEDG(IBL,IS) - USAV(IBL,IS)
      DDS2 = D2_U2*DUE2
C
      CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWAKI,UEI)
      CALL BLKIN
C
C---- check for transition and set TRAN, XT, etc. if found
      IF(TRAN) THEN
        CALL TRCHEK
C       DP mod: check for infinite loop condition
        IF (XFOIL_FAIL) RETURN
        AMI = AMPL2
      ENDIF
C     DP mod: added SILENT_MODE option
      IF(IBL.EQ.ITRAN(IS) .AND. .NOT.TRAN .AND. .NOT.SILENT_MODE) THEN
       WRITE(*,*) 'SETBL: Xtr???  n1 n2: ', AMPL1, AMPL2
      ENDIF
C
C---- assemble 10x4 linearized system for dCtau, dTh, dDs, dUe, dXi
C     at the previous "1" station and the current "2" station
C
      IF(IBL.EQ.IBLTE(IS)+1) THEN
C
C----- define quantities at start of wake, adding TE base thickness to Dstar
       TTE = THET(IBLTE(1),1) + THET(IBLTE(2),2)
       DTE = DSTR(IBLTE(1),1) + DSTR(IBLTE(2),2) + ANTE
       CTE = ( CTAU(IBLTE(1),1)*THET(IBLTE(1),1)
     &       + CTAU(IBLTE(2),2)*THET(IBLTE(2),2) ) / TTE
       CALL TESYS(CTE,TTE,DTE)
C
       TTE_TTE1 = 1.0
       TTE_TTE2 = 1.0
       DTE_MTE1 =               1.0 / UEDG(IBLTE(1),1)
       DTE_UTE1 = -DSTR(IBLTE(1),1) / UEDG(IBLTE(1),1)
       DTE_MTE2 =               1.0 / UEDG(IBLTE(2),2)
       DTE_UTE2 = -DSTR(IBLTE(2),2) / UEDG(IBLTE(2),2)
       CTE_CTE1 = THET(IBLTE(1),1)/TTE
       CTE_CTE2 = THET(IBLTE(2),2)/TTE
       CTE_TTE1 = (CTAU(IBLTE(1),1) - CTE)/TTE
       CTE_TTE2 = (CTAU(IBLTE(2),2) - CTE)/TTE
C
C----- re-define D1 sensitivities wrt m since D1 depends on both TE Ds values
       DO 35 JS=1, 2
         DO 350 JBL=2, NBL(JS)
           J  = IPAN(JBL,JS)
           JV = ISYS(JBL,JS)
           D1_M(JV) = DTE_UTE1*UTE1_M(JV) + DTE_UTE2*UTE2_M(JV)
  350    CONTINUE
   35  CONTINUE
       D1_M(JVTE1) = D1_M(JVTE1) + DTE_MTE1
       D1_M(JVTE2) = D1_M(JVTE2) + DTE_MTE2
C
C----- "forced" changes from  UEDG --- USAV=UINV+dij*MASS  mismatch
       DUE1 = 0.
       DDS1 = DTE_UTE1*(UEDG(IBLTE(1),1) - USAV(IBLTE(1),1))
     &      + DTE_UTE2*(UEDG(IBLTE(2),2) - USAV(IBLTE(2),2))
C
      ELSE
C
       CALL BLSYS
C
      ENDIF
C
C
C---- Save wall shear and equil. max shear coefficient for plotting output
      TAU(IBL,IS) = 0.5*R2*U2*U2*CF2
      DIS(IBL,IS) =     R2*U2*U2*U2*DI2*HS2*0.5
      CTQ(IBL,IS) = CQ2
      DELT(IBL,IS) = DE2
      USLP(IBL,IS) = 1.60/(1.0+US2)
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c      IF(WAKE) THEN
c        ALD = DLCON
c      ELSE
c       ALD = 1.0
c      ENDIF
cC
c      IF(TURB .AND. .NOT.WAKE) THEN
c        GCC = GCCON
c        HKC     = HK2 - 1.0 - GCC/RT2
c        IF(HKC .LT. 0.01) THEN
c         HKC = 0.01
c        ENDIF
c       ELSE
c        HKC = HK2 - 1.0
c       ENDIF
cC
c       HR = HKC     / (GACON*ALD*HK2)
c       UQ = (0.5*CF2 - HR**2) / (GBCON*D2)
cC
c       IF(TURB) THEN
c        IBLP = MIN(IBL+1,NBL(IS))
c        IBLM = MAX(IBL-1,2      )
c        DXSSI = XSSI(IBLP,IS) - XSSI(IBLM,IS)
c        IF(DXXSI.EQ.0.0) DXSSI = 1.0
c        GUXD(IBL,IS) = -LOG(UEDG(IBLP,IS)/UEDG(IBLM,IS)) / DXSSI
c        GUXQ(IBL,IS) = -UQ
c       ELSE
c        GUXD(IBL,IS) = 0.0
c        GUXQ(IBL,IS) = 0.0
c       ENDIF
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C---- set XI sensitivities wrt LE Ue changes
      IF(IS.EQ.1) THEN
       XI_ULE1 =  SST_GO
       XI_ULE2 = -SST_GP
      ELSE
       XI_ULE1 = -SST_GO
       XI_ULE2 =  SST_GP
      ENDIF
C
C---- stuff BL system coefficients into main Jacobian matrix
C
      DO 40 JV=1, NSYS
        VM(1,JV,IV) = VS1(1,3)*D1_M(JV) + VS1(1,4)*U1_M(JV)
     &              + VS2(1,3)*D2_M(JV) + VS2(1,4)*U2_M(JV)
     &              + (VS1(1,5) + VS2(1,5) + VSX(1))
     &               *(XI_ULE1*ULE1_M(JV) + XI_ULE2*ULE2_M(JV))
   40 CONTINUE
C
      VB(1,1,IV) = VS1(1,1)
      VB(1,2,IV) = VS1(1,2)
C
      VA(1,1,IV) = VS2(1,1)
      VA(1,2,IV) = VS2(1,2)
C
      IF(LALFA) THEN
       VDEL(1,2,IV) = VSR(1)*RE_CLMR + VSM(1)*MSQ_CLMR
      ELSE
       VDEL(1,2,IV) = 
     &       (VS1(1,4)*U1_A + VS1(1,3)*D1_A)
     &     + (VS2(1,4)*U2_A + VS2(1,3)*D2_A)
     &     + (VS1(1,5) + VS2(1,5) + VSX(1))
     &      *(XI_ULE1*ULE1_A + XI_ULE2*ULE2_A)
      ENDIF
C
      VDEL(1,1,IV) = VSREZ(1)
     &   + (VS1(1,4)*DUE1 + VS1(1,3)*DDS1)
     &   + (VS2(1,4)*DUE2 + VS2(1,3)*DDS2)
     &   + (VS1(1,5) + VS2(1,5) + VSX(1))
     &    *(XI_ULE1*DULE1 + XI_ULE2*DULE2)
C
C
      DO 50 JV=1, NSYS
        VM(2,JV,IV) = VS1(2,3)*D1_M(JV) + VS1(2,4)*U1_M(JV)
     &              + VS2(2,3)*D2_M(JV) + VS2(2,4)*U2_M(JV)
     &              + (VS1(2,5) + VS2(2,5) + VSX(2))
     &               *(XI_ULE1*ULE1_M(JV) + XI_ULE2*ULE2_M(JV))
   50 CONTINUE
C
      VB(2,1,IV)  = VS1(2,1)
      VB(2,2,IV)  = VS1(2,2)
C
      VA(2,1,IV) = VS2(2,1)
      VA(2,2,IV) = VS2(2,2)
C
      IF(LALFA) THEN
       VDEL(2,2,IV) = VSR(2)*RE_CLMR + VSM(2)*MSQ_CLMR
      ELSE
       VDEL(2,2,IV) = 
     &       (VS1(2,4)*U1_A + VS1(2,3)*D1_A)
     &     + (VS2(2,4)*U2_A + VS2(2,3)*D2_A)
     &     + (VS1(2,5) + VS2(2,5) + VSX(2))
     &      *(XI_ULE1*ULE1_A + XI_ULE2*ULE2_A)
      ENDIF
C
      VDEL(2,1,IV) = VSREZ(2)
     &   + (VS1(2,4)*DUE1 + VS1(2,3)*DDS1)
     &   + (VS2(2,4)*DUE2 + VS2(2,3)*DDS2)
     &   + (VS1(2,5) + VS2(2,5) + VSX(2))
     &    *(XI_ULE1*DULE1 + XI_ULE2*DULE2)
C
C
      DO 60 JV=1, NSYS
        VM(3,JV,IV) = VS1(3,3)*D1_M(JV) + VS1(3,4)*U1_M(JV)
     &              + VS2(3,3)*D2_M(JV) + VS2(3,4)*U2_M(JV)
     &              + (VS1(3,5) + VS2(3,5) + VSX(3))
     &               *(XI_ULE1*ULE1_M(JV) + XI_ULE2*ULE2_M(JV))
   60 CONTINUE
C
      VB(3,1,IV) = VS1(3,1)
      VB(3,2,IV) = VS1(3,2)
C
      VA(3,1,IV) = VS2(3,1)
      VA(3,2,IV) = VS2(3,2)
C
      IF(LALFA) THEN
       VDEL(3,2,IV) = VSR(3)*RE_CLMR + VSM(3)*MSQ_CLMR
      ELSE
       VDEL(3,2,IV) = 
     &       (VS1(3,4)*U1_A + VS1(3,3)*D1_A)
     &     + (VS2(3,4)*U2_A + VS2(3,3)*D2_A)
     &     + (VS1(3,5) + VS2(3,5) + VSX(3))
     &      *(XI_ULE1*ULE1_A + XI_ULE2*ULE2_A)
      ENDIF
C
      VDEL(3,1,IV) = VSREZ(3)
     &   + (VS1(3,4)*DUE1 + VS1(3,3)*DDS1)
     &   + (VS2(3,4)*DUE2 + VS2(3,3)*DDS2)
     &   + (VS1(3,5) + VS2(3,5) + VSX(3))
     &    *(XI_ULE1*DULE1 + XI_ULE2*DULE2)
C
C
      IF(IBL.EQ.IBLTE(IS)+1) THEN
C
C----- redefine coefficients for TTE, DTE, etc
       VZ(1,1)    = VS1(1,1)*CTE_CTE1
       VZ(1,2)    = VS1(1,1)*CTE_TTE1 + VS1(1,2)*TTE_TTE1
       VB(1,1,IV) = VS1(1,1)*CTE_CTE2
       VB(1,2,IV) = VS1(1,1)*CTE_TTE2 + VS1(1,2)*TTE_TTE2
C
       VZ(2,1)    = VS1(2,1)*CTE_CTE1
       VZ(2,2)    = VS1(2,1)*CTE_TTE1 + VS1(2,2)*TTE_TTE1
       VB(2,1,IV) = VS1(2,1)*CTE_CTE2
       VB(2,2,IV) = VS1(2,1)*CTE_TTE2 + VS1(2,2)*TTE_TTE2
C
       VZ(3,1)    = VS1(3,1)*CTE_CTE1
       VZ(3,2)    = VS1(3,1)*CTE_TTE1 + VS1(3,2)*TTE_TTE1
       VB(3,1,IV) = VS1(3,1)*CTE_CTE2
       VB(3,2,IV) = VS1(3,1)*CTE_TTE2 + VS1(3,2)*TTE_TTE2
C
      ENDIF
C
C---- turbulent intervals will follow if currently at transition interval
      IF(TRAN) THEN
        TURB = .TRUE.
C
C------ save transition location
        ITRAN(IS) = IBL
        TFORCE(IS) = TRFORC
        XSSITR(IS) = XT
C
C------ interpolate airfoil geometry to find transition x/c
C-      (for user output)
        IF(IS.EQ.1) THEN
         STR = SST - XT
        ELSE
         STR = SST + XT
        ENDIF
        CHX = XTE - XLE
        CHY = YTE - YLE
        CHSQ = CHX**2 + CHY**2
        XTR = SEVAL(STR,X,XP,S,N)
        YTR = SEVAL(STR,Y,YP,S,N)
        XOCTR(IS) = ((XTR-XLE)*CHX + (YTR-YLE)*CHY)/CHSQ
        YOCTR(IS) = ((YTR-YLE)*CHX - (XTR-XLE)*CHY)/CHSQ
      ENDIF
C
      TRAN = .FALSE.
C
      IF(IBL.EQ.IBLTE(IS)) THEN
C----- set "2" variables at TE to wake correlations for next station
C
       TURB = .TRUE.
       WAKE = .TRUE.
       CALL BLVAR(3)
       CALL BLMID(3)
      ENDIF
C
      DO 80 JS=1, 2
        DO 810 JBL=2, NBL(JS)
          JV = ISYS(JBL,JS)
          U1_M(JV) = U2_M(JV)
          D1_M(JV) = D2_M(JV)
  810   CONTINUE
   80 CONTINUE
C
      U1_A = U2_A
      D1_A = D2_A
C
      DUE1 = DUE2
      DDS1 = DDS2
C      
C---- set BL variables for next station
C     DP mod: to remove need for EQUIVALENCE COM1, COM2, and COMMON
      call com2_to_com1()
C      DO 190 ICOM=1, NCOM
C        COM1(ICOM) = COM2(ICOM)
C  190 CONTINUE
C
C---- next streamwise station
 1000 CONTINUE
C 
C     DP mod: added SILENT_MODE option
      IF (.NOT. SILENT_MODE) THEN
        IF(TFORCE(IS)) THEN
         WRITE(*,9100) IS,XOCTR(IS),ITRAN(IS)
 9100    FORMAT(1X,'Side',I2,' forced transition at x/c = ',F7.4,I5)
        ELSE
         WRITE(*,9200) IS,XOCTR(IS),ITRAN(IS)
 9200    FORMAT(1X,'Side',I2,'  free  transition at x/c = ',F7.4,I5)
        ENDIF
      ENDIF
C
C---- next airfoil side
 2000 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C      Custom solver for coupled viscous-inviscid Newton system:
C
C        A  |  |  .  |  |  .  |    d       R       S
C        B  A  |  .  |  |  .  |    d       R       S
C        |  B  A  .  |  |  .  |    d       R       S
C        .  .  .  .  |  |  .  |    d   =   R - dRe S
C        |  |  |  B  A  |  .  |    d       R       S
C        |  Z  |  |  B  A  .  |    d       R       S
C        .  .  .  .  .  .  .  |    d       R       S
C        |  |  |  |  |  |  B  A    d       R       S
C
C       A, B, Z  3x3  blocks containing linearized BL equation coefficients
C       |        3x1  vectors containing mass defect influence 
C                     coefficients on Ue
C       d        3x1  unknown vectors (Newton deltas for Ctau, Theta, m)
C       R        3x1  residual vectors
C       S        3x1  Re influence vectors
C
C===================================================================70
      SUBROUTINE BLSOLV

      use my_equivalence, only : my_equiv_3_2, my_equiv_3_1
      use xfoil_inc
C
      IVTE1 = ISYS(IBLTE(1),1)
C
      VACC1 = VACCEL
      VACC2 = VACCEL * 2.0 / (S(N) - S(1))
      VACC3 = VACCEL * 2.0 / (S(N) - S(1))
C
      DO 1000 IV=1, NSYS
C
        IVP = IV + 1
C
C====== Invert VA(IV) block
C
C------ normalize first row
        PIVOT = 1.0 / VA(1,1,IV)
        VA(1,2,IV) = VA(1,2,IV) * PIVOT
        DO 10 L=IV, NSYS
          VM(1,L,IV) = VM(1,L,IV)*PIVOT
   10   CONTINUE
        VDEL(1,1,IV) = VDEL(1,1,IV)*PIVOT
        VDEL(1,2,IV) = VDEL(1,2,IV)*PIVOT
C
C------ eliminate lower first column in VA block
        DO 15 K=2, 3
          VTMP = VA(K,1,IV)
          VA(K,2,IV) = VA(K,2,IV) - VTMP*VA(1,2,IV)
          DO 150 L=IV, NSYS
            VM(K,L,IV) = VM(K,L,IV) - VTMP*VM(1,L,IV)
  150     CONTINUE
          VDEL(K,1,IV) = VDEL(K,1,IV) - VTMP*VDEL(1,1,IV)
          VDEL(K,2,IV) = VDEL(K,2,IV) - VTMP*VDEL(1,2,IV)
   15   CONTINUE
C
C
C------ normalize second row
        PIVOT = 1.0 / VA(2,2,IV)
        DO 20 L=IV, NSYS
          VM(2,L,IV) = VM(2,L,IV)*PIVOT
   20   CONTINUE
        VDEL(2,1,IV) = VDEL(2,1,IV)*PIVOT
        VDEL(2,2,IV) = VDEL(2,2,IV)*PIVOT
C
C------ eliminate lower second column in VA block
        K = 3
        VTMP = VA(K,2,IV)
        DO 250 L=IV, NSYS
          VM(K,L,IV) = VM(K,L,IV) - VTMP*VM(2,L,IV)
  250   CONTINUE
        VDEL(K,1,IV) = VDEL(K,1,IV) - VTMP*VDEL(2,1,IV)
        VDEL(K,2,IV) = VDEL(K,2,IV) - VTMP*VDEL(2,2,IV)
C
C
C------ normalize third row
        PIVOT = 1.0/VM(3,IV,IV)
        DO 350 L=IVP, NSYS
          VM(3,L,IV) = VM(3,L,IV)*PIVOT
  350   CONTINUE
        VDEL(3,1,IV) = VDEL(3,1,IV)*PIVOT
        VDEL(3,2,IV) = VDEL(3,2,IV)*PIVOT
C
C
C------ eliminate upper third column in VA block
        VTMP1 = VM(1,IV,IV)
        VTMP2 = VM(2,IV,IV)
        DO 450 L=IVP, NSYS
          VM(1,L,IV) = VM(1,L,IV) - VTMP1*VM(3,L,IV)
          VM(2,L,IV) = VM(2,L,IV) - VTMP2*VM(3,L,IV)
  450   CONTINUE
        VDEL(1,1,IV) = VDEL(1,1,IV) - VTMP1*VDEL(3,1,IV)
        VDEL(2,1,IV) = VDEL(2,1,IV) - VTMP2*VDEL(3,1,IV)
        VDEL(1,2,IV) = VDEL(1,2,IV) - VTMP1*VDEL(3,2,IV)
        VDEL(2,2,IV) = VDEL(2,2,IV) - VTMP2*VDEL(3,2,IV)
C
C------ eliminate upper second column in VA block
        VTMP = VA(1,2,IV)
        DO 460 L=IVP, NSYS
          VM(1,L,IV) = VM(1,L,IV) - VTMP*VM(2,L,IV)
  460   CONTINUE
        VDEL(1,1,IV) = VDEL(1,1,IV) - VTMP*VDEL(2,1,IV)
        VDEL(1,2,IV) = VDEL(1,2,IV) - VTMP*VDEL(2,2,IV)
C
C
        IF(IV.EQ.NSYS) GO TO 1000
C
C====== Eliminate VB(IV+1) block, rows  1 -> 3
        DO 50 K=1, 3
          VTMP1 = VB(K, 1,IVP)
          VTMP2 = VB(K, 2,IVP)
          VTMP3 = VM(K,IV,IVP)
          DO 510 L=IVP, NSYS
            VM(K,L,IVP) = VM(K,L,IVP)
     &        - (  VTMP1*VM(1,L,IV)
     &           + VTMP2*VM(2,L,IV)
     &           + VTMP3*VM(3,L,IV) )
  510     CONTINUE
          VDEL(K,1,IVP) = VDEL(K,1,IVP)
     &        - (  VTMP1*VDEL(1,1,IV)
     &           + VTMP2*VDEL(2,1,IV)
     &           + VTMP3*VDEL(3,1,IV) )
          VDEL(K,2,IVP) = VDEL(K,2,IVP)
     &        - (  VTMP1*VDEL(1,2,IV)
     &           + VTMP2*VDEL(2,2,IV)
     &           + VTMP3*VDEL(3,2,IV) )
   50   CONTINUE
C
        IF(IV.EQ.IVTE1) THEN
C------- eliminate VZ block
         IVZ = ISYS(IBLTE(2)+1,2)
C
         DO 55 K=1, 3
           VTMP1 = VZ(K,1)
           VTMP2 = VZ(K,2)
           DO 515 L=IVP, NSYS
             VM(K,L,IVZ) = VM(K,L,IVZ)
     &         - (  VTMP1*VM(1,L,IV)
     &            + VTMP2*VM(2,L,IV) )
  515      CONTINUE
           VDEL(K,1,IVZ) = VDEL(K,1,IVZ)
     &         - (  VTMP1*VDEL(1,1,IV)
     &            + VTMP2*VDEL(2,1,IV) )
           VDEL(K,2,IVZ) = VDEL(K,2,IVZ)
     &         - (  VTMP1*VDEL(1,2,IV)
     &            + VTMP2*VDEL(2,2,IV) )
   55    CONTINUE
        ENDIF
C
        IF(IVP.EQ.NSYS) GO TO 1000
C
C====== Eliminate lower VM column
        DO 60 KV=IV+2, NSYS
          VTMP1 = VM(1,IV,KV)
          VTMP2 = VM(2,IV,KV)
          VTMP3 = VM(3,IV,KV)
C
          IF(ABS(VTMP1).GT.VACC1) THEN
          DO 610 L=IVP, NSYS
            VM(1,L,KV) = VM(1,L,KV) - VTMP1*VM(3,L,IV)
  610     CONTINUE
          VDEL(1,1,KV) = VDEL(1,1,KV) - VTMP1*VDEL(3,1,IV)
          VDEL(1,2,KV) = VDEL(1,2,KV) - VTMP1*VDEL(3,2,IV)
          ENDIF
C
          IF(ABS(VTMP2).GT.VACC2) THEN
          DO 620 L=IVP, NSYS
            VM(2,L,KV) = VM(2,L,KV) - VTMP2*VM(3,L,IV)
  620     CONTINUE
          VDEL(2,1,KV) = VDEL(2,1,KV) - VTMP2*VDEL(3,1,IV)
          VDEL(2,2,KV) = VDEL(2,2,KV) - VTMP2*VDEL(3,2,IV)
          ENDIF
C
          IF(ABS(VTMP3).GT.VACC3) THEN
          DO 630 L=IVP, NSYS
            VM(3,L,KV) = VM(3,L,KV) - VTMP3*VM(3,L,IV)
  630     CONTINUE
          VDEL(3,1,KV) = VDEL(3,1,KV) - VTMP3*VDEL(3,1,IV)
          VDEL(3,2,KV) = VDEL(3,2,KV) - VTMP3*VDEL(3,2,IV)
          ENDIF
C
   60   CONTINUE
C
 1000 CONTINUE
C
C
C
      DO 2000 IV=NSYS, 2, -1
C
C------ eliminate upper VM columns
        VTMP = VDEL(3,1,IV)
        DO 81 KV=IV-1, 1, -1
          VDEL(1,1,KV) = VDEL(1,1,KV) - VM(1,IV,KV)*VTMP
          VDEL(2,1,KV) = VDEL(2,1,KV) - VM(2,IV,KV)*VTMP
          VDEL(3,1,KV) = VDEL(3,1,KV) - VM(3,IV,KV)*VTMP
   81   CONTINUE
C
        VTMP = VDEL(3,2,IV)
        DO 82 KV=IV-1, 1, -1
          VDEL(1,2,KV) = VDEL(1,2,KV) - VM(1,IV,KV)*VTMP
          VDEL(2,2,KV) = VDEL(2,2,KV) - VM(2,IV,KV)*VTMP
          VDEL(3,2,KV) = VDEL(3,2,KV) - VM(3,IV,KV)*VTMP
   82   CONTINUE
C
 2000 CONTINUE
C
C     DP mod: copy from VA to UNEW to remove need for EQUIVALENCE
C       and from VA to U_AC to remove need for EQUIVALENCe
      DO 83 IV = 1, IVX
        call my_equiv_3_2(VA, UNEW, (/ 1, 1, 1 /), (/ 1, 1 /),
     &                   (/ IV, 1 /), 2)
        call my_equiv_3_2(VA, UNEW, (/ 1, 1, 1 /), (/ 1, 1 /),
     &                   (/ IV, 2 /), 2)
        call my_equiv_3_2(VA, U_AC, (/ 1, 1, IVX /), (/ 1, 1 /),
     &                   (/ IV, 1 /), 2)
        call my_equiv_3_2(VA, U_AC, (/ 1, 1, IVX /), (/ 1, 1 /),
     &                   (/ IV, 2 /), 2)
   83 CONTINUE
C
C     DP mod: copy from VB to QNEW to remove need for EQUIVALENCE
C       and from VB to Q_AC
      DO 84 IV = 1, IQX
        call my_equiv_3_1(VB, QNEW, (/ 1, 1, 1 /), 1, IV, 2)
        call my_equiv_3_1(VB, Q_AC, (/ 1, 1, IVX /), 1, IV, 2)
C
   84 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C      Adds on Newton deltas to boundary layer variables.
C      Checks for excessive changes and underrelaxes if necessary.
C      Calculates max and rms changes.
C      Also calculates the change in the global variable "AC".
C        If LALFA=.TRUE. , "AC" is CL
C        If LALFA=.FALSE., "AC" is alpha
C
C===================================================================70
      SUBROUTINE UPDATE

      use xfoil_inc
C      REAL*8 :: UNEW(IVX,2), U_AC(IVX,2)
C      REAL*8 :: QNEW(IQX),   Q_AC(IQX)
C      EQUIVALENCE (VA(1,1,1), UNEW(1,1)) ,
C     &            (VB(1,1,1), QNEW(1)  )
C      EQUIVALENCE (VA(1,1,IVX), U_AC(1,1)) ,
C     &            (VB(1,1,IVX), Q_AC(1)  )
      REAL*8 :: MSQ
C
C---- max allowable alpha changes per iteration
      DALMAX =  0.5*DTOR
      DALMIN = -0.5*DTOR
C
C---- max allowable CL change per iteration
      DCLMAX =  0.5
      DCLMIN = -0.5
      IF(MATYP.NE.1) DCLMIN = MAX(-0.5 , -0.9*CL)
C
      HSTINV = GAMM1*(MINF/QINF)**2 / (1.0 + 0.5*GAMM1*MINF**2)
C
C---- calculate new Ue distribution assuming no under-relaxation
C-    also set the sensitivity of Ue wrt to alpha or Re
      DO 1 IS=1, 2
        DO 10 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)
C
          DUI    = 0.
          DUI_AC = 0.
          DO 100 JS=1, 2
            DO 1000 JBL=2, NBL(JS)
              J  = IPAN(JBL,JS)
              JV = ISYS(JBL,JS)
              UE_M = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
              DUI    = DUI    + UE_M*(MASS(JBL,JS)+VDEL(3,1,JV))
              DUI_AC = DUI_AC + UE_M*(            -VDEL(3,2,JV))
 1000       CONTINUE
  100     CONTINUE
C
C-------- UINV depends on "AC" only if "AC" is alpha
          IF(LALFA) THEN
           UINV_AC = 0.
          ELSE
           UINV_AC = UINV_A(IBL,IS)
          ENDIF
C
          UNEW(IBL,IS) = UINV(IBL,IS) + DUI
          U_AC(IBL,IS) = UINV_AC      + DUI_AC
C
   10   CONTINUE
    1 CONTINUE
C
C---- set new Qtan from new Ue with appropriate sign change
      DO 2 IS=1, 2
        DO 20 IBL=2, IBLTE(IS)
          I = IPAN(IBL,IS)
          QNEW(I) = VTI(IBL,IS)*UNEW(IBL,IS)
          Q_AC(I) = VTI(IBL,IS)*U_AC(IBL,IS)
   20   CONTINUE
    2 CONTINUE
C
C---- calculate new CL from this new Qtan
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
      CLNEW = 0.
      CL_A  = 0.
      CL_MS = 0.
      CL_AC = 0.
C
      I = 1
      CGINC = 1.0 - (QNEW(I)/QINF)**2
      CPG1  = CGINC / (BETA + BFAC*CGINC)
      CPG1_MS = -CPG1/(BETA + BFAC*CGINC)*(BETA_MSQ + BFAC_MSQ*CGINC)
C
      CPI_Q = -2.0*QNEW(I)/QINF**2
      CPC_CPI = (1.0 - BFAC*CPG1)/ (BETA + BFAC*CGINC)
      CPG1_AC = CPC_CPI*CPI_Q*Q_AC(I)
C
      DO 3 I=1, N
        IP = I+1
        IF(I.EQ.N) IP = 1
C
        CGINC = 1.0 - (QNEW(IP)/QINF)**2
        CPG2  = CGINC / (BETA + BFAC*CGINC)
        CPG2_MS = -CPG2/(BETA + BFAC*CGINC)*(BETA_MSQ + BFAC_MSQ*CGINC)
C
        CPI_Q = -2.0*QNEW(IP)/QINF**2
        CPC_CPI = (1.0 - BFAC*CPG2)/ (BETA + BFAC*CGINC)
        CPG2_AC = CPC_CPI*CPI_Q*Q_AC(IP)
C
        DX   =  (X(IP) - X(I))*CA + (Y(IP) - Y(I))*SA
        DX_A = -(X(IP) - X(I))*SA + (Y(IP) - Y(I))*CA
C
        AG    = 0.5*(CPG2    + CPG1   )
        AG_MS = 0.5*(CPG2_MS + CPG1_MS)
        AG_AC = 0.5*(CPG2_AC + CPG1_AC)
C
        CLNEW = CLNEW + DX  *AG
        CL_A  = CL_A  + DX_A*AG
        CL_MS = CL_MS + DX  *AG_MS
        CL_AC = CL_AC + DX  *AG_AC
C
        CPG1    = CPG2
        CPG1_MS = CPG2_MS
        CPG1_AC = CPG2_AC
    3 CONTINUE
C
C---- initialize under-relaxation factor
      RLX = 1.0
C
      IF(LALFA) THEN
C===== alpha is prescribed: AC is CL
C
C----- set change in Re to account for CL changing, since Re = Re(CL)
       DAC = (CLNEW - CL) / (1.0 - CL_AC - CL_MS*2.0*MINF*MINF_CL)
C
C----- set under-relaxation factor if Re change is too large
       IF(RLX*DAC .GT. DCLMAX) RLX = DCLMAX/DAC
       IF(RLX*DAC .LT. DCLMIN) RLX = DCLMIN/DAC
C
      ELSE
C===== CL is prescribed: AC is alpha
C
C----- set change in alpha to drive CL to prescribed value
       DAC = (CLNEW - CLSPEC) / (0.0 - CL_AC - CL_A)
C
C----- set under-relaxation factor if alpha change is too large
       IF(RLX*DAC .GT. DALMAX) RLX = DALMAX/DAC
       IF(RLX*DAC .LT. DALMIN) RLX = DALMIN/DAC
C
      ENDIF
C
      RMSBL = 0.
      RMXBL = 0.
C
      DHI = 1.5
      DLO = -.5
C
C---- calculate changes in BL variables and under-relaxation if needed
      DO 4 IS=1, 2
        DO 40 IBL=2, NBL(IS)
          IV = ISYS(IBL,IS)
C


C-------- set changes without underrelaxation
          DCTAU = VDEL(1,1,IV) - DAC*VDEL(1,2,IV)
          DTHET = VDEL(2,1,IV) - DAC*VDEL(2,2,IV)
          DMASS = VDEL(3,1,IV) - DAC*VDEL(3,2,IV)
          DUEDG = UNEW(IBL,IS) + DAC*U_AC(IBL,IS)  -  UEDG(IBL,IS)
          DDSTR = (DMASS - DSTR(IBL,IS)*DUEDG)/UEDG(IBL,IS)
C
C-------- normalize changes
          IF(IBL.LT.ITRAN(IS)) DN1 = DCTAU / 10.0
          IF(IBL.GE.ITRAN(IS)) DN1 = DCTAU / CTAU(IBL,IS)
          DN2 = DTHET / THET(IBL,IS)
          DN3 = DDSTR / DSTR(IBL,IS)
          DN4 = ABS(DUEDG)/0.25
C
C-------- accumulate for rms change
          RMSBL = RMSBL + DN1**2 + DN2**2 + DN3**2 + DN4**2
C          
C-------- see if Ctau needs underrelaxation
          RDN1 = RLX*DN1
          IF(ABS(DN1) .GT. ABS(RMXBL)) THEN
           RMXBL = DN1
           IF(IBL.LT.ITRAN(IS)) VMXBL = 'n'
           IF(IBL.GE.ITRAN(IS)) VMXBL = 'C'
           IMXBL = IBL
           ISMXBL = IS
          ENDIF
          IF(RDN1 .GT. DHI) RLX = DHI/DN1
          IF(RDN1 .LT. DLO) RLX = DLO/DN1
C
C-------- see if Theta needs underrelaxation
          RDN2 = RLX*DN2
          IF(ABS(DN2) .GT. ABS(RMXBL)) THEN
           RMXBL = DN2
           VMXBL = 'T'
           IMXBL = IBL
           ISMXBL = IS
          ENDIF
          IF(RDN2 .GT. DHI) RLX = DHI/DN2
          IF(RDN2 .LT. DLO) RLX = DLO/DN2
C
C-------- see if Dstar needs underrelaxation
          RDN3 = RLX*DN3
          IF(ABS(DN3) .GT. ABS(RMXBL)) THEN
           RMXBL = DN3
           VMXBL = 'D'
           IMXBL = IBL
           ISMXBL = IS
          ENDIF
          IF(RDN3 .GT. DHI) RLX = DHI/DN3
          IF(RDN3 .LT. DLO) RLX = DLO/DN3
C
C-------- see if Ue needs underrelaxation
          RDN4 = RLX*DN4
          IF(ABS(DN4) .GT. ABS(RMXBL)) THEN
           RMXBL = DUEDG
           VMXBL = 'U'
           IMXBL = IBL
           ISMXBL = IS
          ENDIF
          IF(RDN4 .GT. DHI) RLX = DHI/DN4
          IF(RDN4 .LT. DLO) RLX = DLO/DN4
C
   40   CONTINUE
    4 CONTINUE
C
C---- set true rms change
      RMSBL = SQRT( RMSBL / (4.0*FLOAT( NBL(1)+NBL(2) )) )
C
C
      IF(LALFA) THEN
C----- set underrelaxed change in Reynolds number from change in lift
       CL = CL + RLX*DAC
      ELSE
C----- set underrelaxed change in alpha
       ALFA = ALFA + RLX*DAC
       ADEG = ALFA/DTOR
      ENDIF
C
C---- update BL variables with underrelaxed changes
      DO 5 IS=1, 2
        DO 50 IBL=2, NBL(IS)
          IV = ISYS(IBL,IS)
C
          DCTAU = VDEL(1,1,IV) - DAC*VDEL(1,2,IV)
          DTHET = VDEL(2,1,IV) - DAC*VDEL(2,2,IV)
          DMASS = VDEL(3,1,IV) - DAC*VDEL(3,2,IV)
          DUEDG = UNEW(IBL,IS) + DAC*U_AC(IBL,IS)  -  UEDG(IBL,IS)
          DDSTR = (DMASS - DSTR(IBL,IS)*DUEDG)/UEDG(IBL,IS)
C
          CTAU(IBL,IS) = CTAU(IBL,IS) + RLX*DCTAU
          THET(IBL,IS) = THET(IBL,IS) + RLX*DTHET
          DSTR(IBL,IS) = DSTR(IBL,IS) + RLX*DDSTR
          UEDG(IBL,IS) = UEDG(IBL,IS) + RLX*DUEDG
C
          IF(IBL.GT.IBLTE(IS)) THEN
           IW = IBL - IBLTE(IS)
           DSWAKI = WGAP(IW)
          ELSE
           DSWAKI = 0.
          ENDIF
C
C-------- eliminate absurd transients
          IF(IBL.GE.ITRAN(IS))
     &      CTAU(IBL,IS) = MIN( CTAU(IBL,IS) , 0.25 )
C
          IF(IBL.LE.IBLTE(IS)) THEN
            HKLIM = 1.02
          ELSE
            HKLIM = 1.00005
          ENDIF
          MSQ = UEDG(IBL,IS)**2*HSTINV
     &        / (GAMM1*(1.0 - 0.5*UEDG(IBL,IS)**2*HSTINV))
          DSW = DSTR(IBL,IS) - DSWAKI
          CALL DSLIM(DSW,THET(IBL,IS),UEDG(IBL,IS),MSQ,HKLIM)
          DSTR(IBL,IS) = DSW + DSWAKI
C
C-------- set new mass defect (nonlinear update)
          MASS(IBL,IS) = DSTR(IBL,IS) * UEDG(IBL,IS)
C
   50   CONTINUE
C
C------ make sure there are no "islands" of negative Ue
        DO IBL = 3, IBLTE(IS)
          IF(UEDG(IBL-1,IS) .GT. 0.0 .AND.
     &       UEDG(IBL  ,IS) .LE. 0.0       ) THEN
           UEDG(IBL,IS) = UEDG(IBL-1,IS)
           MASS(IBL,IS) = DSTR(IBL,IS) * UEDG(IBL,IS)
          ENDIF
        ENDDO
    5 CONTINUE
C
C
C---- equate upper wake arrays to lower wake arrays
      DO 6 KBL=1, NBL(2)-IBLTE(2)
        CTAU(IBLTE(1)+KBL,1) = CTAU(IBLTE(2)+KBL,2)
        THET(IBLTE(1)+KBL,1) = THET(IBLTE(2)+KBL,2)
        DSTR(IBLTE(1)+KBL,1) = DSTR(IBLTE(2)+KBL,2)
        UEDG(IBLTE(1)+KBL,1) = UEDG(IBLTE(2)+KBL,2)
         TAU(IBLTE(1)+KBL,1) =  TAU(IBLTE(2)+KBL,2)
         DIS(IBLTE(1)+KBL,1) =  DIS(IBLTE(2)+KBL,2)
         CTQ(IBLTE(1)+KBL,1) =  CTQ(IBLTE(2)+KBL,2)
        DELT(IBLTE(1)+KBL,1) = DELT(IBLTE(2)+KBL,2)
        TSTR(IBLTE(1)+KBL,1) = TSTR(IBLTE(2)+KBL,2)
    6 CONTINUE
C
      RETURN
      END

C===================================================================70
C
C      Copies "2" variables to "1" variables to remove need for
C      EQUIVALENCE and COMMON blocks
C
C===================================================================70
      subroutine com2_to_com1()

      use xbl_inc
      implicit none

      X1     = X2    
      U1     = U2    
      T1     = T2    
      D1     = D2    
      S1     = S2    
      AMPL1  = AMPL2 
      U1_UEI = U2_UEI 
      U1_MS  = U2_MS  
      DW1    = DW2   
      H1     = H2    
      H1_T1  = H2_T2 
      H1_D1  = H2_D2 
      M1     = M2    
      M1_U1  = M2_U2 
      M1_MS  = M2_MS 
      R1     = R2    
      R1_U1  = R2_U2 
      R1_MS  = R2_MS 
      V1     = V2    
      V1_U1  = V2_U2  
      V1_MS  = V2_MS  
      V1_RE  = V2_RE  
      HK1    = HK2    
      HK1_U1 = HK2_U2  
      HK1_T1 = HK2_T2  
      HK1_D1 = HK2_D2  
      HK1_MS = HK2_MS  
      HS1    = HS2     
      HS1_U1 = HS2_U2    
      HS1_T1 = HS2_T2   
      HS1_D1 = HS2_D2  
      HS1_MS = HS2_MS  
      HS1_RE = HS2_RE  
      HC1    = HC2    
      HC1_U1 = HC2_U2   
      HC1_T1 = HC2_T2   
      HC1_D1 = HC2_D2   
      HC1_MS = HC2_MS   
      RT1    = RT2    
      RT1_U1 = RT2_U2  
      RT1_T1 = RT2_T2  
      RT1_MS = RT2_MS  
      RT1_RE = RT2_RE  
      CF1    = CF2    
      CF1_U1 = CF2_U2  
      CF1_T1 = CF2_T2  
      CF1_D1 = CF2_D2  
      CF1_MS = CF2_MS  
      CF1_RE = CF2_RE  
      DI1    = DI2    
      DI1_U1 = DI2_U2  
      DI1_T1 = DI2_T2  
      DI1_D1 = DI2_D2  
      DI1_S1 = DI2_S2  
      DI1_MS = DI2_MS  
      DI1_RE = DI2_RE  
      US1    = US2    
      US1_U1 = US2_U2  
      US1_T1 = US2_T2  
      US1_D1 = US2_D2  
      US1_MS = US2_MS  
      US1_RE = US2_RE  
      CQ1    = CQ2    
      CQ1_U1 = CQ2_U2  
      CQ1_T1 = CQ2_T2  
      CQ1_D1 = CQ2_D2  
      CQ1_MS = CQ2_MS  
      CQ1_RE = CQ2_RE  
      DE1    = DE2    
      DE1_U1 = DE2_U2  
      DE1_T1 = DE2_T2  
      DE1_D1 = DE2_D2  
      DE1_MS = DE2_MS  

      end subroutine com2_to_com1

C===================================================================70
C
C      Copies "1" variables to C1SAV array
C
C===================================================================70
      subroutine store_c1sav()

      use xbl_inc
      implicit none

      C1SAV(1)  = X1
      C1SAV(2)  = U1
      C1SAV(3)  = T1
      C1SAV(4)  = D1
      C1SAV(5)  = S1
      C1SAV(6)  = AMPL1
      C1SAV(7)  = U1_UEI
      C1SAV(8)  = U1_MS
      C1SAV(9)  = DW1
      C1SAV(10) = H1
      C1SAV(11) = H1_T1
      C1SAV(12) = H1_D1
      C1SAV(13) = M1
      C1SAV(14) = M1_U1
      C1SAV(15) = M1_MS
      C1SAV(16) = R1
      C1SAV(17) = R1_U1
      C1SAV(18) = R1_MS
      C1SAV(19) = V1
      C1SAV(20) = V1_U1
      C1SAV(21) = V1_MS
      C1SAV(22) = V1_RE
      C1SAV(23) = HK1
      C1SAV(24) = HK1_U1
      C1SAV(25) = HK1_T1
      C1SAV(26) = HK1_D1
      C1SAV(27) = HK1_MS
      C1SAV(28) = HS1
      C1SAV(29) = HS1_U1
      C1SAV(30) = HS1_T1
      C1SAV(31) = HS1_D1
      C1SAV(32) = HS1_MS
      C1SAV(33) = HS1_RE
      C1SAV(34) = HC1
      C1SAV(35) = HC1_U1
      C1SAV(36) = HC1_T1
      C1SAV(37) = HC1_D1
      C1SAV(38) = HC1_MS
      C1SAV(39) = RT1
      C1SAV(40) = RT1_U1
      C1SAV(41) = RT1_T1
      C1SAV(42) = RT1_MS
      C1SAV(43) = RT1_RE
      C1SAV(44) = CF1
      C1SAV(45) = CF1_U1
      C1SAV(46) = CF1_T1
      C1SAV(47) = CF1_D1
      C1SAV(48) = CF1_MS
      C1SAV(49) = CF1_RE
      C1SAV(50) = DI1
      C1SAV(51) = DI1_U1
      C1SAV(52) = DI1_T1
      C1SAV(53) = DI1_D1
      C1SAV(54) = DI1_S1
      C1SAV(55) = DI1_MS
      C1SAV(56) = DI1_RE
      C1SAV(57) = US1
      C1SAV(58) = US1_U1
      C1SAV(59) = US1_T1
      C1SAV(60) = US1_D1
      C1SAV(61) = US1_MS
      C1SAV(62) = US1_RE
      C1SAV(63) = CQ1
      C1SAV(64) = CQ1_U1
      C1SAV(65) = CQ1_T1
      C1SAV(66) = CQ1_D1
      C1SAV(67) = CQ1_MS
      C1SAV(68) = CQ1_RE
      C1SAV(69) = DE1
      C1SAV(70) = DE1_U1
      C1SAV(71) = DE1_T1
      C1SAV(72) = DE1_D1
      C1SAV(73) = DE1_MS

      end subroutine store_c1sav

C===================================================================70
C
C      Copies "1" variables from C1SAV array
C
C===================================================================70
      subroutine from_c1sav()

      use xbl_inc
      implicit none

      X1     = C1SAV(1)  
      U1     = C1SAV(2)    
      T1     = C1SAV(3)    
      D1     = C1SAV(4)      
      S1     = C1SAV(5)    
      AMPL1  = C1SAV(6)   
      U1_UEI = C1SAV(7)  
      U1_MS  = C1SAV(8)  
      DW1    = C1SAV(9)  
      H1     = C1SAV(10)   
      H1_T1  = C1SAV(11)  
      H1_D1  = C1SAV(12)  
      M1     = C1SAV(13)  
      M1_U1  = C1SAV(14)   
      M1_MS  = C1SAV(15) 
      R1     = C1SAV(16)  
      R1_U1  = C1SAV(17)
      R1_MS  = C1SAV(18)
      V1     = C1SAV(19)   
      V1_U1  = C1SAV(20) 
      V1_MS  = C1SAV(21)    
      V1_RE  = C1SAV(22)    
      HK1    = C1SAV(23)   
      HK1_U1 = C1SAV(24)    
      HK1_T1 = C1SAV(25)      
      HK1_D1 = C1SAV(26)        
      HK1_MS = C1SAV(27)       
      HS1    = C1SAV(28)    
      HS1_U1 = C1SAV(29)     
      HS1_T1 = C1SAV(30)     
      HS1_D1 = C1SAV(31)   
      HS1_MS = C1SAV(32)     
      HS1_RE = C1SAV(33)    
      HC1    = C1SAV(34)  
      HC1_U1 = C1SAV(35)   
      HC1_T1 = C1SAV(36)     
      HC1_D1 = C1SAV(37) 
      HC1_MS = C1SAV(38) 
      RT1    = C1SAV(39)  
      RT1_U1 = C1SAV(40) 
      RT1_T1 = C1SAV(41) 
      RT1_MS = C1SAV(42) 
      RT1_RE = C1SAV(43) 
      CF1    = C1SAV(44) 
      CF1_U1 = C1SAV(45) 
      CF1_T1 = C1SAV(46) 
      CF1_D1 = C1SAV(47) 
      CF1_MS = C1SAV(48) 
      CF1_RE = C1SAV(49) 
      DI1    = C1SAV(50)   
      DI1_U1 = C1SAV(51) 
      DI1_T1 = C1SAV(52) 
      DI1_D1 = C1SAV(53) 
      DI1_S1 = C1SAV(54) 
      DI1_MS = C1SAV(55) 
      DI1_RE = C1SAV(56) 
      US1    = C1SAV(57)   
      US1_U1 = C1SAV(58) 
      US1_T1 = C1SAV(59) 
      US1_D1 = C1SAV(60) 
      US1_MS = C1SAV(61) 
      US1_RE = C1SAV(62) 
      CQ1    = C1SAV(63)  
      CQ1_U1 = C1SAV(64) 
      CQ1_T1 = C1SAV(65) 
      CQ1_D1 = C1SAV(66) 
      CQ1_MS = C1SAV(67) 
      CQ1_RE = C1SAV(68) 
      DE1    = C1SAV(69)  
      DE1_U1 = C1SAV(70) 
      DE1_T1 = C1SAV(71) 
      DE1_D1 = C1SAV(72) 
      DE1_MS = C1SAV(73) 

      end subroutine from_c1sav

C===================================================================70
C
C      Copies "2" variables to C2SAV array
C
C===================================================================70
      subroutine store_c2sav()

      use xbl_inc
      implicit none

      C2SAV(1)  = X2
      C2SAV(2)  = U2
      C2SAV(3)  = T2
      C2SAV(4)  = D2
      C2SAV(5)  = S2
      C2SAV(6)  = AMPL2
      C2SAV(7)  = U2_UEI
      C2SAV(8)  = U2_MS
      C2SAV(9)  = DW2
      C2SAV(10) = H2
      C2SAV(11) = H2_T2
      C2SAV(12) = H2_D2
      C2SAV(13) = M2
      C2SAV(14) = M2_U2
      C2SAV(15) = M2_MS
      C2SAV(16) = R2
      C2SAV(17) = R2_U2
      C2SAV(18) = R2_MS
      C2SAV(19) = V2
      C2SAV(20) = V2_U2
      C2SAV(21) = V2_MS
      C2SAV(22) = V2_RE
      C2SAV(23) = HK2
      C2SAV(24) = HK2_U2
      C2SAV(25) = HK2_T2
      C2SAV(26) = HK2_D2
      C2SAV(27) = HK2_MS
      C2SAV(28) = HS2
      C2SAV(29) = HS2_U2
      C2SAV(30) = HS2_T2
      C2SAV(31) = HS2_D2
      C2SAV(32) = HS2_MS
      C2SAV(33) = HS2_RE
      C2SAV(34) = HC2
      C2SAV(35) = HC2_U2
      C2SAV(36) = HC2_T2
      C2SAV(37) = HC2_D2
      C2SAV(38) = HC2_MS
      C2SAV(39) = RT2
      C2SAV(40) = RT2_U2
      C2SAV(41) = RT2_T2
      C2SAV(42) = RT2_MS
      C2SAV(43) = RT2_RE
      C2SAV(44) = CF2
      C2SAV(45) = CF2_U2
      C2SAV(46) = CF2_T2
      C2SAV(47) = CF2_D2
      C2SAV(48) = CF2_MS
      C2SAV(49) = CF2_RE
      C2SAV(50) = DI2
      C2SAV(51) = DI2_U2
      C2SAV(52) = DI2_T2
      C2SAV(53) = DI2_D2
      C2SAV(54) = DI2_S2
      C2SAV(55) = DI2_MS
      C2SAV(56) = DI2_RE
      C2SAV(57) = US2
      C2SAV(58) = US2_U2
      C2SAV(59) = US2_T2
      C2SAV(60) = US2_D2
      C2SAV(61) = US2_MS
      C2SAV(62) = US2_RE
      C2SAV(63) = CQ2
      C2SAV(64) = CQ2_U2
      C2SAV(65) = CQ2_T2
      C2SAV(66) = CQ2_D2
      C2SAV(67) = CQ2_MS
      C2SAV(68) = CQ2_RE
      C2SAV(69) = DE2
      C2SAV(70) = DE2_U2
      C2SAV(71) = DE2_T2
      C2SAV(72) = DE2_D2
      C2SAV(73) = DE2_MS

      end subroutine store_c2sav

C===================================================================70
C
C      Copies "2" variables from C2SAV array
C
C===================================================================70
      subroutine from_c2sav()

      use xbl_inc
      implicit none

      X2     = C2SAV(1) 
      U2     = C2SAV(2) 
      T2     = C2SAV(3) 
      D2     = C2SAV(4) 
      S2     = C2SAV(5) 
      AMPL2  = C2SAV(6) 
      U2_UEI = C2SAV(7) 
      U2_MS  = C2SAV(8) 
      DW2    = C2SAV(9) 
      H2     = C2SAV(10)
      H2_T2  = C2SAV(11)
      H2_D2  = C2SAV(12)
      M2     = C2SAV(13)
      M2_U2  = C2SAV(14)
      M2_MS  = C2SAV(15)
      R2     = C2SAV(16)
      R2_U2  = C2SAV(17)
      R2_MS  = C2SAV(18)
      V2     = C2SAV(19)
      V2_U2  = C2SAV(20)
      V2_MS  = C2SAV(21)
      V2_RE  = C2SAV(22)
      HK2    = C2SAV(23)
      HK2_U2 = C2SAV(24)
      HK2_T2 = C2SAV(25)
      HK2_D2 = C2SAV(26)
      HK2_MS = C2SAV(27)
      HS2    = C2SAV(28)
      HS2_U2 = C2SAV(29)
      HS2_T2 = C2SAV(30)
      HS2_D2 = C2SAV(31)
      HS2_MS = C2SAV(32)
      HS2_RE = C2SAV(33)
      HC2    = C2SAV(34)
      HC2_U2 = C2SAV(35)
      HC2_T2 = C2SAV(36)
      HC2_D2 = C2SAV(37)
      HC2_MS = C2SAV(38)
      RT2    = C2SAV(39)
      RT2_U2 = C2SAV(40)
      RT2_T2 = C2SAV(41)
      RT2_MS = C2SAV(42)
      RT2_RE = C2SAV(43)
      CF2    = C2SAV(44)
      CF2_U2 = C2SAV(45)
      CF2_T2 = C2SAV(46)
      CF2_D2 = C2SAV(47)
      CF2_MS = C2SAV(48)
      CF2_RE = C2SAV(49)
      DI2    = C2SAV(50)
      DI2_U2 = C2SAV(51)
      DI2_T2 = C2SAV(52)
      DI2_D2 = C2SAV(53)
      DI2_S2 = C2SAV(54)
      DI2_MS = C2SAV(55)
      DI2_RE = C2SAV(56)
      US2    = C2SAV(57)
      US2_U2 = C2SAV(58)
      US2_T2 = C2SAV(59)
      US2_D2 = C2SAV(60)
      US2_MS = C2SAV(61)
      US2_RE = C2SAV(62)
      CQ2    = C2SAV(63)
      CQ2_U2 = C2SAV(64)
      CQ2_T2 = C2SAV(65)
      CQ2_D2 = C2SAV(66)
      CQ2_MS = C2SAV(67)
      CQ2_RE = C2SAV(68)
      DE2    = C2SAV(69)
      DE2_U2 = C2SAV(70)
      DE2_T2 = C2SAV(71)
      DE2_D2 = C2SAV(72)
      DE2_MS = C2SAV(73)

      end subroutine from_c2sav
