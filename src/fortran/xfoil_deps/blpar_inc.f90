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

!  Copyright (C) 2017-2019 Daniel Prosser

module blpar_inc

  implicit none
!
!-    SCCON  =  shear coefficient lag constant
!-    GACON  =  G-beta locus constants...
!-    GBCON  =   G = GACON * sqrt(1.0 + GBCON*beta) 
!-    GCCON  =         + GCCON / [H*Rtheta*sqrt(Cf/2)]   <-- wall term
!-    DLCON  =  wall/wake dissipation length ratio  Lo/L
!-    CTCON  =  Ctau weighting coefficient (implied by G-beta constants)
!
  REAL*8 :: SCCON, GACON, GBCON, GCCON, DLCON, CTRCON, CTRCEX, DUXCON, CTCON
  !$omp threadprivate(SCCON, GACON, GBCON, GCCON, DLCON, CTRCON, CTRCEX)
  !$omp threadprivate(DUXCON, CTCON)
  REAL*8 :: CFFAC
  !$omp threadprivate(CFFAC)

end module blpar_inc
