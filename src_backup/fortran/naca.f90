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

!  Copyright (C) 2017 Daniel Prosser

module naca

! Creates naca 4, 5, 6, and 6A series airfoils

  implicit none

  type naca_options_type

    character(3) :: family        ! '4', '4M', '5', '63', '64', '65', '66',
                                  ! '67', '63A', '64A', or '65A'
    double precision :: maxt      ! Max thickness/chord
    double precision :: xmaxt     ! Location of maxt. Only for 4M.
    double precision :: maxc      ! Max camber/chord. Only for 4 and 4M.
    double precision :: xmaxc     ! Location of maxc. All except for 6 and 6A.
    double precision :: design_cl ! Design Cl. Only for 5, 6, and 6A.
    double precision :: a         ! Extent of constant load in x/c. Only for 6.
    double precision :: leidx     ! Leading edge index. Only for 4M.
    logical :: reflexed           ! Whether mean line is reflexed. Only for 5.

  end type naca_options_type

  contains

!=============================================================================80
!
! Cosine point spacing distribution function for a vector with given number of
! points and endpoints.
!
!=============================================================================80
subroutine cosine_spacing(x, nx, x0, xf)

  implicit none

  integer, intent(in) :: nx
  double precision, dimension(nx), intent(inout) :: x
  double precision, intent(in) :: x0, xf

  double precision :: r, xm, pi, theta
  integer :: i

  pi = acos(-1.0)
  r = 0.5*(xf - x0)
  xm = 0.5*(x0 + xf)
  do i = 1, nx
    theta = pi * (1.0 - real(i-1)/real(nx-1))
    x(i) = xm + r*cos(theta)
  end do

end subroutine cosine_spacing

!=============================================================================80
!
! Creates 4, 5, 6, or 6A series seed airfoil
!
!=============================================================================80
subroutine naca_456(naca_options, pointsmcl, foil)

  use vardef,          only : airfoil_type
  use memory_util,     only : allocate_airfoil
  use NacaAuxiliary

  type(naca_options_type), intent(in) :: naca_options
  integer, intent(in) :: pointsmcl
  type(airfoil_type), intent(out) :: foil

  integer :: i
  character(3) :: family
  double precision, dimension(pointsmcl) :: x, zm, zmp, zt, xu, xl, zu, zl

! Get mean camber line x-coordinates via normal distribution spacing

  call cosine_spacing(x, pointsmcl, 0.d0, 1.d0)

! Thickness profile

  family = adjustl(naca_options%family)
  select case (adjustl(family))
    case('4', '5')
      call Thickness4(naca_options%maxt, x, zt)
    case('4M', '4m')
      call Thickness4M(naca_options%maxt, naca_options%leidx,                  &
                       naca_options%xmaxt, x, zt)
    case('63')
      call Thickness6(1, naca_options%maxt, x, zt)
    case('64')
      call Thickness6(2, naca_options%maxt, x, zt)
    case('65')
      call Thickness6(3, naca_options%maxt, x, zt)
    case('66')
      call Thickness6(4, naca_options%maxt, x, zt)
    case('67')
      call Thickness6(5, naca_options%maxt, x, zt)
    case('63A', '63a')
      call Thickness6(6, naca_options%maxt, x, zt)
    case('64A', '64a')
      call Thickness6(7, naca_options%maxt, x, zt)
    case('65A', '65a')
      call Thickness6(8, naca_options%maxt, x, zt)
    case default
      write(*,*) family//" is not a valid family."
      stop
  end select

! Mean camber line

  select case (adjustl(family))
    case ('4', '4M', '4m')
      call MeanLine2(naca_options%maxc, naca_options%xmaxc, x, zm, zmp)
    case ('5')
      if (naca_options%reflexed) then
        call MeanLine3Reflex(naca_options%design_cl, naca_options%xmaxc, x, zm,&
                             zmp)
      else
        call MeanLine3(naca_options%design_cl, naca_options%xmaxc, x, zm, zmp)
      end if
    case ('63', '64', '65', '66', '67')
      call MeanLine6(naca_options%a, naca_options%design_cl, x, zm, zmp)
    case ('63A', '63a', '64A', '64a', '65A', '65a')
      call MeanLine6M(naca_options%design_cl, x, zm, zmp)
    case default
      write(*,*) family//" is not a valid family."
      stop
  end select

! Combine thickness and camber

  call CombineThicknessAndCamber(x, zt, zm, zmp, xu, zu, xl, zl)

! Format coordinates in a single loop (in airfoil_type derived type)

  foil%npoint = 2*pointsmcl - 1
  call allocate_airfoil(foil)
  do i = 1, pointsmcl
    foil%x(i) = xu(pointsmcl-i+1)
    foil%z(i) = zu(pointsmcl-i+1)
  end do
  do i = 1, pointsmcl - 1
    foil%x(i+pointsmcl) = xl(i+1)
    foil%z(i+pointsmcl) = zl(i+1)
  end do

end subroutine naca_456

end module naca
