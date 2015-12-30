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

module my_equivalence

! Provides routines to copy array values so that equivalence statements can
! be removed.  This is needed for parallelization.

  implicit none

  integer, parameter :: toA = 1
  integer, parameter :: toB = 2

  contains

!=============================================================================80
!
! Manual equivalence between A (3-D) and B (2-D)
! Copies B(bcopyidx(1),bcopyidx(2)) to the appropriate location in A, if 
!   direction == toA (1).  Alternatively, copies from the appropriate location
!   in A to B(bcopyidx(1),bcopyidx(2)), if direction == toB (2).
!
! Used to replace an equivalence statement like the following:
!   EQUIVALENCE(A(aequividx(1),aequividx(2),aequividx(3)),
!               B(bequividx(1),bequividx(2))) 
!
!=============================================================================80
subroutine my_equiv_3_2(A, B, aequividx, bequividx, bcopyidx, direction)

  double precision, dimension(:,:,:), intent(inout) :: A
  double precision, dimension(:,:), intent(inout) :: B
  integer, dimension(3), intent(in) :: aequividx
  integer, dimension(2), intent(in) :: bequividx, bcopyidx
  integer, intent(in) :: direction

  integer na1, na2, na3, nb1, nb2
  integer an, bn
  integer aequiv, bequiv
  integer, dimension(3) :: aidx

  na1 = size(A,1)
  na2 = size(A,2)
  na3 = size(A,3)
  nb1 = size(B,1)
  nb2 = size(B,2)

! Find virtual 1-D equivalence points

  aequiv = equiv_pt_3d(na1, na2, aequividx)
  bequiv = equiv_pt_2d(nb1, bequividx)

! 1-D transformation of point to copy from B

  bn = equiv_pt_2d(nb1, bcopyidx)

! 1-D point in A to copy to

  an = bn - bequiv + aequiv

! Translation to 3-D index in A

  aidx = inv_equiv_pt_3d(an, na1, na2)

! Copy data

  if (direction == toA) then
    A(aidx(1),aidx(2),aidx(3)) = B(bcopyidx(1),bcopyidx(2))
  else
    B(bcopyidx(1),bcopyidx(2)) = A(aidx(1),aidx(2),aidx(3))
  end if

end subroutine my_equiv_3_2

!=============================================================================80
!
! Manual equivalence of B (2-D) and A (3-D)
! Copies A(acopyidx(1),acopyidx(2),acopyidx(3)) to the appropriate location in 
!   B, if direction == toB (2).  Alternatively, copies to 
!   A(acopyidx(1),acopyidx(2),acopyidx(3)) from the appropriate location in B,
!   if direction == toA (1).
!
! Used to replace an equivalence statement like the following:
!   EQUIVALENCE(A(aequividx(1),aequividx(2),aequividx(3)),
!               B(bequividx(1),bequividx(2))) 
!
!=============================================================================80
subroutine my_equiv_2_3(A, B, aequividx, bequividx, acopyidx, direction)

  double precision, dimension(:,:,:), intent(inout) :: A
  double precision, dimension(:,:), intent(inout) :: B
  integer, dimension(3), intent(in) :: aequividx, acopyidx
  integer, dimension(2), intent(in) :: bequividx
  integer, intent(in) :: direction

  integer na1, na2, na3, nb1, nb2
  integer an, bn
  integer aequiv, bequiv
  integer, dimension(2) :: bidx

  na1 = size(A,1)
  na2 = size(A,2)
  na3 = size(A,3)
  nb1 = size(B,1)
  nb2 = size(B,2)

! Find virtual 1-D equivalence points

  aequiv = equiv_pt_3d(na1, na2, aequividx)
  bequiv = equiv_pt_2d(nb1, bequividx)

! 1-D transformation of point to copy from a

  an = equiv_pt_3d(na1, na2, acopyidx)

! 1-D point in B to copy to

  bn = an - aequiv + bequiv

! Translation to 2-D index in B

  bidx = inv_equiv_pt_2d(bn, nb1)

! Copy data

  if (direction == toB) then
    B(bidx(1),bidx(2)) = A(acopyidx(1),acopyidx(2),acopyidx(3))
  else
    A(acopyidx(1),acopyidx(2),acopyidx(3)) = B(bidx(1),bidx(2))
  end if

end subroutine my_equiv_2_3

!=============================================================================80
!
! Manual equivalence between A (3-D) and B (1-D)
! Copies B(bcopyidx) to the appropriate location in A, if 
!   direction == toA (1).  Alternatively, copies from the appropriate location
!   in A to B(bcopyidx), if direction == toB (2).
!
! Used to replace an equivalence statement like the following:
!   EQUIVALENCE(A(aequividx(1),aequividx(2),aequividx(3)),
!               B(bequividx) 
!
!=============================================================================80
subroutine my_equiv_3_1(A, B, aequividx, bequividx, bcopyidx, direction)

  double precision, dimension(:,:,:), intent(inout) :: A
  double precision, dimension(:), intent(inout) :: B
  integer, dimension(3), intent(in) :: aequividx
  integer, intent(in) :: bequividx, bcopyidx, direction

  integer na1, na2, na3, nb1
  integer an, bn
  integer aequiv, bequiv
  integer, dimension(3) :: aidx

  na1 = size(A,1)
  na2 = size(A,2)
  na3 = size(A,3)
  nb1 = size(B,1)

! Find virtual 1-D equivalence points

  aequiv = equiv_pt_3d(na1, na2, aequividx)
  bequiv = bequividx

! 1-D transformation of point to copy from B

  bn = bcopyidx

! 1-D point in A to copy to

  an = bn - bequiv + aequiv

! Translation to 3-D index in A

  aidx = inv_equiv_pt_3d(an, na1, na2)

! Copy data

  if (direction == toA) then
    A(aidx(1),aidx(2),aidx(3)) = B(bcopyidx)
  else
    B(bcopyidx) = A(aidx(1),aidx(2),aidx(3))
  end if

end subroutine my_equiv_3_1

!=============================================================================80
!
! Determines 1-D point in memory for a 3D array with dimensions n1, n2, n3
! for the index (aidx(1), aidx(2), aidx(3))
!
!=============================================================================80
function equiv_pt_3d(n1, n2, aidx) result(equiv_pt)

  integer, intent(in) :: n1, n2
  integer, dimension(3), intent(in) :: aidx
  integer :: equiv_pt

  equiv_pt = (aidx(3)-1)*n1*n2 + (aidx(2)-1)*n1 + aidx(1)

end function equiv_pt_3d

!=============================================================================80
!
! Determines 3-D indices corresponding to 1-D point in memory (for 3-D array)
!
!=============================================================================80
function inv_equiv_pt_3d(pt, n1, n2) result(inv_equiv_pt)

  integer, intent(in) :: pt, n1, n2
  integer, dimension(3) :: inv_equiv_pt

  inv_equiv_pt(3) = ceiling(dble(pt)/dble(n1*n2))
  inv_equiv_pt(2) = ceiling(dble(pt - (inv_equiv_pt(3)-1)*n1*n2)/dble(n1))
  inv_equiv_pt(1) = pt - (inv_equiv_pt(2)-1)*n1 - (inv_equiv_pt(3)-1)*n1*n2

end function inv_equiv_pt_3d

!=============================================================================80
!
! Determines 1-D point in memory for a 2D array with dimensions n1, n2 for the
! index (aidx(1), aidx(2))
!
!=============================================================================80
function equiv_pt_2d(n1, aidx) result(equiv_pt)

  integer, intent(in) :: n1
  integer, dimension(2), intent(in) :: aidx
  integer :: equiv_pt

  equiv_pt = (aidx(2)-1)*n1 + aidx(1)

end function equiv_pt_2d

!=============================================================================80
!
! Determines 2-D indices corresponding to 1-D point in memory (for 2-D array)
!
!=============================================================================80
function inv_equiv_pt_2d(pt, n1) result(inv_equiv_pt)

  integer, intent(in) :: pt, n1
  integer, dimension(2) :: inv_equiv_pt

  inv_equiv_pt(2) = ceiling(dble(pt)/dble(n1))
  inv_equiv_pt(1) = pt - (inv_equiv_pt(2)-1)*n1

end function inv_equiv_pt_2d

end module my_equivalence
