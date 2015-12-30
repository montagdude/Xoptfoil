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

module math_deps

! Contains various math functions and numerical methods

  implicit none

  contains

!=============================================================================80
!
! Function to get x = inv(A)*C using gaussian elimination
!
!=============================================================================80
function lmult(A,C) result(X)

  double precision, dimension(:,:), intent(in) :: A
  double precision, dimension(:), intent(in) :: C
  double precision, dimension(size(C,1)) :: X
  double precision, dimension(size(C,1),size(C,1)+1) :: Q
  integer :: N, i, j, R
  double precision :: elim, pivot, rscale, rsum, eps
  eps = 1D-16

! Initialize

  N = size(C,1)
  if (size(A,1) /= N .or. size(A,2) /= N) then
    write(*,*)
    write(*,*) 'Error: for A*X = C and size(C) = Nx1, size(A) must be NxN'
    write(*,*)
    stop
  end if
  X(:) =  0.d0
  Q(:,1:N) = A(:,:)
  Q(:,N+1) = C(:)

! Gaussian elimination loop to put in upper triangular form

  do R = 1, N-1
    pivot = Q(R,R)
    do i = R+1, N
      elim = Q(i,R)
      if (abs(elim) > eps) then
        rscale = elim/pivot
        Q(i,:) = Q(i,:) - rscale*Q(R,:)
      end if
    end do
  end do

! Solution loop

  do i = N, 1, -1
    rsum = Q(i,N+1)
    do j = N, i+1, -1
      if (abs(Q(i,j)) > eps) rsum = rsum - Q(i,j)*X(j)
    end do
    if (Q(i,i) == 0) then
      write(*,*)
      write(*,*) 'Error in lmult: singular matrix.'
      stop
    else
      X(i) = rsum/Q(i,i)
    end if
  end do

end function lmult

!=============================================================================80
!
! Normal distribution function, used for small spacing at ends and greater in
! the middle
!
!=============================================================================80
function normal_dist(x, sig, mu) result(val)

  double precision, intent(in) :: x, sig, mu
  double precision val, pi

  pi = acos(-1.d0)
  val = 1.d0/(sig*sqrt(2.d0*pi))*exp(-(x-mu)**2.d0/(2.d0*sig**2.d0))

end function normal_dist

!=============================================================================80
!
! Vector norm (since not all compilers may include it by default)
!
!=============================================================================80
function norm_2(vector) result(val)

  double precision, dimension(:), intent(in) :: vector
  double precision :: val
  integer :: nelem, i

! Determine size

  nelem = size(vector)

! Get vector norm

  val = 0.d0
  do i = 1, nelem
    val = val + vector(i)**2.d0
  end do
  val = sqrt(val)

end function norm_2

!=============================================================================80
!
! Interpolates a vector y with original coordinates x to a new set of
! coordinates xnew
!
!=============================================================================80
subroutine interp_vector(x, y, xnew, ynew)

  double precision, dimension(:), intent(in) :: x, y, xnew
  double precision, dimension(:), intent(inout) :: ynew

  logical :: isbtwn
  integer :: i, pt1, npt, nptnew

  npt = size(x,1)
  nptnew = size(xnew,1)

  pt1 = 1
  do i = 1, nptnew

!   Find interpolants

    isbtwn = .false.
    do while (.not. isbtwn .and. (pt1 < npt))
      isbtwn = between(x(pt1), xnew(i), x(pt1+1))
      if (.not. isbtwn) then
        pt1 = pt1 + 1
        if (pt1 == npt) then
          write(*,*)
          write(*,*) 'Warning: could not find interpolants.'
          write(*,*) 'x: ', xnew(i), 'xmax: ', x(npt)
          stop
        end if
      end if
    end do

!   Interpolate points

    ynew(i) = interp1(x(pt1), x(pt1+1), xnew(i), y(pt1), y(pt1+1))

  end do

end subroutine interp_vector

!=============================================================================80
!
! Interpolates between two points
!
!=============================================================================80
function interp1(x1, x2, x, y1, y2) result(y)

  double precision, intent(in) :: x1, x2, x, y1, y2
  double precision y

  y = y1 + (y2 - y1)*(x - x1)/(x2 - x1)

end function interp1

!=============================================================================80
!
! Determines if B is between A and C
!
!=============================================================================80
function between(A, B, C) result(test)

  double precision, intent(in) :: A, B, C
  logical test

  if ((B >= A) .and. (B <= C)) then
    test = .true.
  else
    test = .false.
  end if 

end function between

!=============================================================================80
!
! Computes curvature for a function gam(s) = x(s) + y(s)
!
!=============================================================================80
function curvature(npt, x, y)

  integer, intent(in) :: npt
  double precision, dimension(npt), intent(in) :: x, y
  double precision, dimension(npt) :: curvature

  integer :: i
  double precision, dimension(npt) :: svec
  double precision :: se, se2
  double precision :: xe, ye, xe2, ye2
  double precision :: xs, ys, xs2, ys2

! Airfoil length vector s 

  svec(1) = 0.d0
  do i = 2, npt
    svec(i) = svec(i-1) + sqrt((x(i)-x(i-1))**2.d0 + (y(i)-y(i-1))**2.d0)
  end do

! Compute first and second derivatives and curvature vector

  do i = 1, npt

    if (i == 1) then

!     Grid metric ds/de and d2s/de2

      se = derv1f(svec(i+2), svec(i+1), svec(i), 1.d0)
      se2 = derv2f(svec(i+2), svec(i+1), svec(i), 1.d0)

!     Derivatives of x and y with respect to the grid parameter e

      xe = derv1f(x(i+2), x(i+1), x(i), 1.d0)
      ye = derv1f(y(i+2), y(i+1), y(i), 1.d0)
      xe2 = derv2f(x(i+2), x(i+1), x(i), 1.d0)
      ye2 = derv2f(y(i+2), y(i+1), y(i), 1.d0)

    elseif (i == npt) then

!     Grid metric ds/de and d2s de2

      se = derv1b(svec(i-2), svec(i-1), svec(i), 1.d0)
      se2 = derv2b(svec(i-2), svec(i-1), svec(i), 1.d0)

!     Derivatives of x and y with respect to the grid parameter e

      xe = derv1b(x(i-2), x(i-1), x(i), 1.d0)
      ye = derv1b(y(i-2), y(i-1), y(i), 1.d0)
      xe2 = derv2b(x(i-2), x(i-1), x(i), 1.d0)
      ye2 = derv2b(y(i-2), y(i-1), y(i), 1.d0)
      
    else

!     Grid metric ds/de and d2s de2

      se = derv1c(svec(i+1), svec(i-1), 1.d0)
      se2 = derv2c(svec(i+1), svec(i), svec(i-1), 1.d0)

!     Derivatives of x and y with respect to the grid parameter e

      xe = derv1c(x(i+1), x(i-1), 1.d0)
      ye = derv1c(y(i+1), y(i-1), 1.d0)
      xe2 = derv2c(x(i+1), x(i), x(i-1), 1.d0)
      ye2 = derv2c(y(i+1), y(i), y(i-1), 1.d0)

    end if

!   Derivatives of x and y with respect to surface length s

    xs = 1.d0/se * xe
    ys = 1.d0/se * ye
    xs2 = 1.d0/se**2.d0 * (xe2 - se2/se*xe)
    ys2 = 1.d0/se**2.d0 * (ye2 - se2/se*ye)

!   Curvature

    curvature(i) = (xs*ys2 - ys*xs2) / (xs**2.d0 + ys**2.d0)**1.5d0

  end do

end function curvature

!=============================================================================80
!
! Forward difference approximation for first derivative (2nd order)
!
!=============================================================================80
function derv1f(u_plus2, u_plus1, u, h)

  double precision, intent(in) :: u_plus2, u_plus1, u, h
  double precision :: derv1f

  derv1f = (-3.d0*u + 4.d0*u_plus1 - u_plus2) / (2.d0*h)

end function derv1f

!=============================================================================80
!
! Backward difference approximation for first derivative (2nd order)
!
!=============================================================================80
function derv1b(u_minus2, u_minus1, u, h)

  double precision, intent(in) :: u_minus2, u_minus1, u, h
  double precision :: derv1b

  derv1b = (3.d0*u - 4.d0*u_minus1 + u_minus2) / (2.d0*h)

end function derv1b

!=============================================================================80
!
! Central difference approximation for first derivative (2nd order)
!
!=============================================================================80
function derv1c(u_plus, u_minus, h)

  double precision, intent(in) :: u_plus, u_minus, h
  double precision :: derv1c

  derv1c = (u_plus - u_minus) / (2.d0*h)

end function derv1c

!=============================================================================80
!
! Forward difference approximation for second-order derivative
!
!=============================================================================80
function derv2f(u_plus2, u_plus, u, h)

  double precision, intent(in) :: u_plus2, u_plus, u, h
  double precision :: derv2f

  derv2f = (u - 2.d0*u_plus + u_plus2) / h**2.d0

end function derv2f

!=============================================================================80
!
! Backward difference approximation for second-order derivative
!
!=============================================================================80
function derv2b(u_minus2, u_minus, u, h)

  double precision, intent(in) :: u_minus2, u_minus, u, h
  double precision :: derv2b

  derv2b = (u - 2.d0*u_minus + u_minus2) / h**2.d0

end function derv2b

!=============================================================================80
!
! Central difference approximation for second-order derivative
!
!=============================================================================80
function derv2c(u_plus, u, u_minus, h)

  double precision, intent(in) :: u_plus, u, u_minus, h
  double precision :: derv2c

  derv2c = (u_plus - 2.d0*u + u_minus) / h**2.d0

end function derv2c

end module math_deps
