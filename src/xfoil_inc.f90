!  This file is part of XOPTFOIL.

!  XOPTFOIL is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.

!  XOPTFOIL is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.

!  You should have received a copy of the GNU General Public License
!  along with XOPTFOIL.  If not, see <http://www.gnu.org/licenses/>.

!  Copyright (C) 2014 -- 2015 Daniel Prosser

module xfoil_inc

! Global include file for xfoil

!------ Primary dimensioning limit parameters
! IQX   number of surface panel nodes + 6
! IWX   number of wake panel nodes
! IPX   number of Qspec(s) distributions
! ISX   number of airfoil sides
!
!------ Derived dimensioning limit parameters
! IBX   number of buffer airfoil nodes
! IMX   number of complex mapping coefficients  Cn
! IZX   number of panel nodes (airfoil + wake)
! IVX   number of nodes along BL on one side of airfoil and wake
! NAX   number of points in stored polar
! NPX   number of polars and reference polars
! NFX   number of points in one reference polar
! NTX   number of points in thickness/camber arrays

  INTEGER, PARAMETER :: IQX=360
  INTEGER, PARAMETER :: ISX=2
  INTEGER, PARAMETER :: IBX=4*IQX
  INTEGER, PARAMETER :: IWX=IQX/8+2
  INTEGER, PARAMETER :: IZX=IQX+IWX
  INTEGER, PARAMETER :: IVX=IQX/2 + IWX + 50
  LOGICAL :: SILENT_MODE
  REAL*8 :: PI, HOPI, QOPI, DTOR
  REAL*8 :: GAM(IQX), GAMU(IQX,2), QINVU(IZX,2), GAM_A(IQX)
  REAL*8 :: QINV(IZX), QINV_A(IZX)
  LOGICAL :: LGAMU, LQAIJ, SHARP, LVISC, LWAKE, LVCONV, LWDIJ, LIPAN, LBLINI
  LOGICAL :: LADIJ, LALFA
  INTEGER :: RETYP, MATYP, ITMAX
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: AIJ, BIJ, DIJ, CIJ
  REAL*8 :: DZDG(IQX), DZDN(IQX), DQDG(IQX), DZDM(IZX), DQDM(IZX)
  REAL*8 :: X(IZX), Y(IZX), NX(IZX), NY(IZX), S(IZX), APANEL(IZX), SIG(IZX)
  REAL*8 :: XP(IZX), YP(IZX)
  INTEGER :: N, NB, NPAN, NW, IST, NSYS
  REAL*8 :: PSIO, QINF, ALFA, Z_QINF, Z_ALFA, Z_QDOF0, Z_QDOF1, Z_QDOF2, Z_QDOF3
  REAL*8 :: ANTE, ASTE, DSTE, ADEG, AMAX
  REAL*8 :: QF0(IQX), QF1(IQX), QF2(IQX), QF3(IQX)
  INTEGER :: AIJPIV(IQX), IBLTE(ISX), NBL(ISX)
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: IPAN, ISYS
  REAL*8 :: SIGTE, GAMTE, SIGTE_A, GAMTE_A, MINF, MINF1, REINF, REINF1
  REAL*8 :: TKLAM, TKL_MSQ, CPSTAR, QSTAR, GAMMA, GAMM1
  REAL*8 :: XCMREF, YCMREF, CL, CM, CD, CDP, CDF, CL_ALF, CL_MSQ, SBLE
  REAL*8 :: XB(IBX), YB(IBX), SB(IBX), XBP(IBX), YBP(IBX), SNEW(5*IBX)
  REAL*8, DIMENSION(:), ALLOCATABLE :: W1, W2, W3, W4, W5, W6
  REAL*8 :: XLE, YLE, XTE, YTE, CHORD, SLE
  REAL*8 :: CVPAR, CTERAT, CTRRAT, XSREF1, XSREF2, XPREF1, XPREF2
  REAL*8 :: MINF_CL, COSA, SINA, ACRIT, RLX, VACCEL
  REAL*8 :: CPI(IZX), CPV(IZX), QVIS(IZX)
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: VTI, XSSI
  REAL*8 :: AWAKE, AVISC, MVISC, CLSPEC, QTAN1, QTAN2, SST, SST_GO, SST_GP
  REAL*8 :: WGAP(IWX), XSTRIP(ISX), XSSITR(ISX) 
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: UINV, UINV_A, UEDG
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: THET, DSTR, CTAU
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: MASS, TAU, DIS, CTQ
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: DELT, TSTR, USLP
  INTEGER :: IDAMP
  INTEGER :: ITRAN(ISX), IMXBL, ISMXBL
  LOGICAL :: TFORCE(ISX)
  REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: VM, VA, VB, VDEL
  REAL*8 :: VZ(3,2), XOCTR(ISX), YOCTR(ISX)
  REAL*8 :: RMSBL, RMXBL, WAKLEN
  REAL*8 :: UNEW(IVX,2), U_AC(IVX,2)
  REAL*8 :: QNEW(IQX), Q_AC(IQX)
  CHARACTER(1) :: VMXBL
  REAL*8 :: THICKB, XTHICKB, THICKM, XTHICKM, CAMBR, XCAMBR 
  LOGICAL :: XFOIL_FAIL

! Make all these variables threadprivate for OpenMP

  !$omp threadprivate(GAM, GAMU, QINVU, GAM_A)
  !$omp threadprivate(QINV, QINV_A)
  !$omp threadprivate(LGAMU, LQAIJ, SHARP, LVISC, LWAKE, LVCONV, LWDIJ, LIPAN)
  !$omp threadprivate(LBLINI, LADIJ, LALFA)
  !$omp threadprivate(RETYP, MATYP, ITMAX)
  !$omp threadprivate(AIJ, BIJ, DIJ, CIJ)
  !$omp threadprivate(DZDG, DZDN, DQDG, DZDM, DQDM)
  !$omp threadprivate(X, Y, NX, NY, S, APANEL, SIG)
  !$omp threadprivate(XP, YP)
  !$omp threadprivate(N, NB, NPAN, NW, IST, NSYS)
  !$omp threadprivate(PSIO, QINF, ALFA, Z_QINF, Z_ALFA, Z_QDOF0, Z_QDOF1)
  !$omp threadprivate(Z_QDOF2, Z_QDOF3)
  !$omp threadprivate(ANTE, ASTE, DSTE, ADEG, AMAX)
  !$omp threadprivate(QF0, QF1, QF2, QF3)
  !$omp threadprivate(AIJPIV, IPAN, IBLTE, NBL, ISYS)
  !$omp threadprivate(SIGTE, GAMTE, SIGTE_A, GAMTE_A, MINF, MINF1)
  !$omp threadprivate(REINF, REINF1)
  !$omp threadprivate(TKLAM, TKL_MSQ, CPSTAR, QSTAR, GAMMA, GAMM1)
  !$omp threadprivate(XCMREF, YCMREF, CL, CM, CD, CDP, CDF)
  !$omp threadprivate(CL_ALF, CL_MSQ, SBLE)
  !$omp threadprivate(XB, YB, SB, XBP, YBP, SNEW)
  !$omp threadprivate(W1, W2, W3, W4, W5, W6)
  !$omp threadprivate(XLE, YLE, XTE, YTE, CHORD, SLE)
  !$omp threadprivate(CVPAR, CTERAT, CTRRAT, XSREF1, XSREF2, XPREF1, XPREF2)
  !$omp threadprivate(MINF_CL, COSA, SINA, ACRIT, RLX, VACCEL)
  !$omp threadprivate(CPI, CPV, QVIS, VTI, XSSI)
  !$omp threadprivate(AWAKE, AVISC, MVISC, CLSPEC, QTAN1, QTAN2)
  !$omp threadprivate(SST, SST_GO, SST_GP)
  !$omp threadprivate(WGAP, UINV, UINV_A, UEDG)
  !$omp threadprivate(THET, DSTR, CTAU, XSTRIP)
  !$omp threadprivate(MASS, TAU, DIS, CTQ)
  !$omp threadprivate(DELT, TSTR, XSSITR, USLP)
  !$omp threadprivate(IDAMP)
  !$omp threadprivate(ITRAN, IMXBL, ISMXBL)
  !$omp threadprivate(TFORCE)
  !$omp threadprivate(VM, VA, VB, VDEL, VZ)
  !$omp threadprivate(XOCTR, YOCTR)
  !$omp threadprivate(RMSBL, RMXBL, WAKLEN)
  !$omp threadprivate(UNEW, U_AC)
  !$omp threadprivate(QNEW, Q_AC)
  !$omp threadprivate(VMXBL)
  !$omp threadprivate(THICKB, XTHICKB, THICKM, XTHICKM, CAMBR, XCAMBR) 
  !$omp threadprivate(XFOIL_FAIL)

end module xfoil_inc
