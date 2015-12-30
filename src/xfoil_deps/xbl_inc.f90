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

module xbl_inc

  implicit none
      
  INTEGER, PARAMETER :: NCOM = 73
  INTEGER :: IDAMPV
  !$omp threadprivate(IDAMPV)
  LOGICAL :: SIMI, TRAN, TURB, WAKE
  !$omp threadprivate(SIMI, TRAN, TURB, WAKE)
  LOGICAL :: TRFORC, TRFREE
  !$omp threadprivate(TRFORC, TRFREE)

  REAL*8 :: X1, U1, T1, D1, S1, AMPL1, U1_UEI, U1_MS, DW1
  !$omp threadprivate(X1, U1, T1, D1, S1, AMPL1, U1_UEI, U1_MS, DW1)
  REAL*8 :: H1, H1_T1, H1_D1
  !$omp threadprivate(H1, H1_T1, H1_D1)
  REAL*8 :: M1, M1_U1, M1_MS
  !$omp threadprivate(M1, M1_U1, M1_MS)
  REAL*8 :: R1, R1_U1, R1_MS
  !$omp threadprivate(R1, R1_U1, R1_MS)
  REAL*8 :: V1, V1_U1, V1_MS, V1_RE
  !$omp threadprivate(V1, V1_U1, V1_MS, V1_RE)
  REAL*8 :: HK1, HK1_U1, HK1_T1, HK1_D1, HK1_MS
  !$omp threadprivate(HK1, HK1_U1, HK1_T1, HK1_D1, HK1_MS)
  REAL*8 :: HS1, HS1_U1, HS1_T1, HS1_D1, HS1_MS, HS1_RE
  !$omp threadprivate(HS1, HS1_U1, HS1_T1, HS1_D1, HS1_MS, HS1_RE)
  REAL*8 :: HC1, HC1_U1, HC1_T1, HC1_D1, HC1_MS
  !$omp threadprivate(HC1, HC1_U1, HC1_T1, HC1_D1, HC1_MS)
  REAL*8 :: RT1, RT1_U1, RT1_T1, RT1_MS, RT1_RE
  !$omp threadprivate(RT1, RT1_U1, RT1_T1, RT1_MS, RT1_RE)
  REAL*8 :: CF1, CF1_U1, CF1_T1, CF1_D1, CF1_MS, CF1_RE
  !$omp threadprivate(CF1, CF1_U1, CF1_T1, CF1_D1, CF1_MS, CF1_RE)
  REAL*8 :: DI1, DI1_U1, DI1_T1, DI1_D1, DI1_S1, DI1_MS, DI1_RE
  !$omp threadprivate(DI1, DI1_U1, DI1_T1, DI1_D1, DI1_S1, DI1_MS, DI1_RE)
  REAL*8 :: US1, US1_U1, US1_T1, US1_D1, US1_MS, US1_RE
  !$omp threadprivate(US1, US1_U1, US1_T1, US1_D1, US1_MS, US1_RE)
  REAL*8 :: CQ1, CQ1_U1, CQ1_T1, CQ1_D1, CQ1_MS, CQ1_RE 
  !$omp threadprivate(CQ1, CQ1_U1, CQ1_T1, CQ1_D1, CQ1_MS, CQ1_RE) 
  REAL*8 :: DE1, DE1_U1, DE1_T1, DE1_D1, DE1_MS            
  !$omp threadprivate(DE1, DE1_U1, DE1_T1, DE1_D1, DE1_MS)
  REAL*8 :: X2, U2, T2, D2, S2, AMPL2, U2_UEI, U2_MS, DW2
  !$omp threadprivate(X2, U2, T2, D2, S2, AMPL2, U2_UEI, U2_MS, DW2)
  REAL*8 :: H2, H2_T2, H2_D2
  !$omp threadprivate(H2, H2_T2, H2_D2)
  REAL*8 :: M2, M2_U2, M2_MS
  !$omp threadprivate(M2, M2_U2, M2_MS)
  REAL*8 :: R2, R2_U2, R2_MS
  !$omp threadprivate(R2, R2_U2, R2_MS)
  REAL*8 :: V2, V2_U2, V2_MS, V2_RE
  !$omp threadprivate(V2, V2_U2, V2_MS, V2_RE)
  REAL*8 :: HK2, HK2_U2, HK2_T2, HK2_D2, HK2_MS
  !$omp threadprivate(HK2, HK2_U2, HK2_T2, HK2_D2, HK2_MS)
  REAL*8 :: HS2, HS2_U2, HS2_T2, HS2_D2, HS2_MS, HS2_RE
  !$omp threadprivate(HS2, HS2_U2, HS2_T2, HS2_D2, HS2_MS, HS2_RE)
  REAL*8 :: HC2, HC2_U2, HC2_T2, HC2_D2, HC2_MS
  !$omp threadprivate(HC2, HC2_U2, HC2_T2, HC2_D2, HC2_MS)
  REAL*8 :: RT2, RT2_U2, RT2_T2, RT2_MS, RT2_RE
  !$omp threadprivate(RT2, RT2_U2, RT2_T2, RT2_MS, RT2_RE)
  REAL*8 :: CF2, CF2_U2, CF2_T2, CF2_D2, CF2_MS, CF2_RE
  !$omp threadprivate(CF2, CF2_U2, CF2_T2, CF2_D2, CF2_MS, CF2_RE)
  REAL*8 :: DI2, DI2_U2, DI2_T2, DI2_D2, DI2_S2, DI2_MS, DI2_RE
  !$omp threadprivate(DI2, DI2_U2, DI2_T2, DI2_D2, DI2_S2, DI2_MS, DI2_RE)
  REAL*8 :: US2, US2_U2, US2_T2, US2_D2, US2_MS, US2_RE
  !$omp threadprivate(US2, US2_U2, US2_T2, US2_D2, US2_MS, US2_RE)
  REAL*8 :: CQ2, CQ2_U2, CQ2_T2, CQ2_D2, CQ2_MS, CQ2_RE
  !$omp threadprivate(CQ2, CQ2_U2, CQ2_T2, CQ2_D2, CQ2_MS, CQ2_RE)
  REAL*8 :: DE2, DE2_U2, DE2_T2, DE2_D2, DE2_MS            
  !$omp threadprivate(DE2, DE2_U2, DE2_T2, DE2_D2, DE2_MS)

  REAL*8 :: CFM, CFM_MS, CFM_RE
  !$omp threadprivate(CFM, CFM_MS, CFM_RE)
  REAL*8 :: CFM_U1, CFM_T1, CFM_D1, CFM_U2, CFM_T2, CFM_D2
  !$omp threadprivate(CFM_U1, CFM_T1, CFM_D1, CFM_U2, CFM_T2, CFM_D2)
  REAL*8 :: XT, XT_A1, XT_MS, XT_RE, XT_XF
  !$omp threadprivate(XT, XT_A1, XT_MS, XT_RE, XT_XF)
  REAL*8 :: XT_X1, XT_T1, XT_D1, XT_U1, XT_X2, XT_T2, XT_D2, XT_U2         
  !$omp threadprivate(XT_X1, XT_T1, XT_D1, XT_U1, XT_X2, XT_T2, XT_D2, XT_U2)

  REAL*8, DIMENSION(NCOM) :: C1SAV, C2SAV
  !$omp threadprivate(C1SAV, C2SAV)
 
  REAL*8 :: DWTE, QINFBL, TKBL, TKBL_MS, RSTBL, RSTBL_MS, HSTINV, HSTINV_MS
  !$omp threadprivate(DWTE, QINFBL, TKBL, TKBL_MS, RSTBL, RSTBL_MS)
  !$omp threadprivate(HSTINV, HSTINV_MS)
  REAL*8 :: REYBL, REYBL_MS, REYBL_RE, GAMBL, GM1BL
  !$omp threadprivate(REYBL, REYBL_MS, REYBL_RE, GAMBL, GM1BL)
  REAL*8 :: HVRA, BULE, XIFORC, AMCRIT          
  !$omp threadprivate(HVRA, BULE, XIFORC, AMCRIT)
 
  REAL*8, DIMENSION(4,5) :: VS1, VS2
  !$omp threadprivate(VS1, VS2)
  REAL*8, DIMENSION(4) :: VSREZ, VSR, VSM, VSX
  !$omp threadprivate(VSREZ, VSR, VSM, VSX)

end module xbl_inc
