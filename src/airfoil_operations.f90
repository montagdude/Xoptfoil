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

!  Copyright (C) 2014 -- 2016 Daniel Prosser

module airfoil_operations

! Performs transformations and other operations on airfoils

  implicit none

! Coefficients for 5th-order polynomial (curve fit for leading edge)

  double precision, dimension(6) :: polynomial_coefs

  contains

!=============================================================================80
!
! Driver subroutine to read or create a seed airfoil
!
!=============================================================================80
subroutine get_seed_airfoil(seed_airfoil, airfoil_file, naca_digits, foil)

  use vardef,       only : airfoil_type
  use xfoil_driver, only : smooth_paneling

  character(*), intent(in) :: seed_airfoil, airfoil_file
  character(4), intent(in) :: naca_digits
  type(airfoil_type), intent(out) :: foil

  type(airfoil_type) :: tempfoil
  integer :: pointsmcl

  if (trim(seed_airfoil) == 'from_file') then

!   Read seed airfoil from file

    call load_airfoil(airfoil_file, tempfoil)

  elseif (trim(seed_airfoil) == 'four_digit') then

!   Create seed airfoil from NACA 4-digit series

    pointsmcl = 200
    call naca_four_digit(naca_digits, pointsmcl, tempfoil)

  else

    write(*,*) "Error: seed_airfoil should be 'from_file' or 'four_digit'."
    write(*,*)
    stop

  end if

! Use Xfoil to smooth airfoil paneling

  call smooth_paneling(tempfoil, 200, foil)

! Calculate leading edge information

  call le_find(foil%x, foil%z, foil%leclose, foil%xle, foil%zle,               &
               foil%addpoint_loc)

! Translate and scale

  call transform_airfoil(foil)

end subroutine get_seed_airfoil

!=============================================================================80
!
! Reads an airfoil from a file, loads it into the airfoil_type, sets ordering
! correctly
!
!=============================================================================80
subroutine load_airfoil(filename, foil)

  use vardef, only : airfoil_type

  character(*), intent(in) :: filename
  type(airfoil_type), intent(out) :: foil

  logical :: labeled

  write(*,*) 'Reading airfoil from file: '//trim(filename)//' ...'
  write(*,*)

! Read number of points and allocate coordinates

  call airfoil_points(filename, foil%npoint, labeled)
  call allocate_airfoil(foil)

! Read airfoil from file

  call airfoil_read(filename, foil%npoint, labeled, foil%x, foil%z)

! Change point ordering to counterclockwise, if necessary

  call cc_ordering(foil)

end subroutine load_airfoil

!=============================================================================80
!
! Creates a NACA 4-digit airfoil (top and bottom surfaces)
!
!=============================================================================80
subroutine naca_four_digit(naca_digits, pointsmcl, foil)

  use parameterization, only : normal_spacing
  use vardef,           only : airfoil_type

  character(4), intent(in) :: naca_digits
  integer, intent(in) :: pointsmcl
  type(airfoil_type), intent(out) :: foil

  integer :: i
  double precision :: chord, d0, width, m, p, t
  double precision, dimension(pointsmcl) :: xmcl, zmcl, zth, th, xt, xb, zt, zb

  write(*,*) 'Generating seed airfoil: NACA '//trim(naca_digits)//' ...'
  write(*,*)

  chord = 1.d0

! Check to make sure input is sane

  do i = 1, 4
    if (.not. isnum(naca_digits(i:i))) &
      call my_stop("NACA digits should be numeric.")
  end do

! Read naca digits into real values

  read(naca_digits(1:1),*) m
  read(naca_digits(2:2),*) p
  read(naca_digits(3:4),*) t
  m = m / 100.d0
  p = p / 10.d0
  t = t / 100.d0

! To prevent divide by zero errors

  if (p == 0) p = 0.3d0

! Make spacing near leading and trailing edge 0.25 x uniform spacing

  d0 = 0.0625*chord/dble(pointsmcl-1)

! Get mean camber line x-coordinates via normal distribution spacing

  xmcl(1) = 0.d0
  do i = 1, pointsmcl - 2
    call normal_spacing(width, chord, i, pointsmcl-1, d0)
    xmcl(i+1) = xmcl(i) + width
  end do
  xmcl(pointsmcl) = 1.d0
  
! Get mean camber line z-coordinates and inclination angle

  do i = 1, pointsmcl
    if (xmcl(i) <= p) then
      zmcl(i) = m/p**2.d0*xmcl(i)*(2.d0*p - xmcl(i))
      th(i) = atan(2.d0*m/p**2.d0*(p - xmcl(i)))
    else
      zmcl(i) = m*(1.d0 - xmcl(i))/(1.d0 - p)**2.d0*(1.d0 + xmcl(i) - 2.d0*p)
      th(i) = atan(2.d0*m/(1.d0 - p)**2.d0*(p - xmcl(i)))
    end if
  end do

! Thickness equation

  do i = 1, pointsmcl
    zth(i) = t/0.2d0*(0.2969d0*sqrt(xmcl(i)) - 0.126d0*xmcl(i) -               &
                      0.3516d0*xmcl(i)**2.d0 + 0.2843d0*xmcl(i)**3.d0 -        &
                      0.1015d0*xmcl(i)**4.d0)
  end do

! Upper and lower surfaces

  do i = 1, pointsmcl
    xt(i) = xmcl(i) - zth(i)*sin(th(i))
    zt(i) = zmcl(i) + zth(i)*cos(th(i))
    xb(i) = xmcl(i) + zth(i)*sin(th(i))
    zb(i) = zmcl(i) - zth(i)*cos(th(i))
  end do

! Format coordinates in a single loop (in airfoil_type derived type)

  foil%npoint = 2*pointsmcl - 1
  call allocate_airfoil(foil)
  do i = 1, pointsmcl
    foil%x(i) = xt(pointsmcl-i+1)
    foil%z(i) = zt(pointsmcl-i+1)
  end do
  do i = 1, pointsmcl - 1
    foil%x(i+pointsmcl) = xb(i+1)
    foil%z(i+pointsmcl) = zb(i+1)
  end do

end subroutine naca_four_digit

!=============================================================================80
!
! Subroutine to get number of points from an airfoil file and to determine
! whether it is labeled or plain.
!
!=============================================================================80
subroutine airfoil_points(filename, npoints, labeled)

  character(*), intent(in) :: filename
  integer, intent(out) :: npoints
  logical, intent(out) :: labeled

  integer :: iunit, ioerr
  double precision :: dummyx, dummyz

! Open airfoil file

  iunit = 12
  open(unit=iunit, file=filename, status='old', position='rewind', iostat=ioerr)
  if (ioerr /= 0) then
     write(*,*) 'Error: cannot find airfoil file '//trim(filename)
     write(*,*)
     stop
  end if

! Read first line; determine if it is a title or not

  read(iunit,*,iostat=ioerr) dummyx, dummyz
  if (ioerr == 0) then
    npoints = 1
    labeled = .false.
  else
    npoints = 0
    labeled = .true.
  end if
  
! Read the rest of the lines

  do 
    read(iunit,*,end=500)
    npoints = npoints + 1
  end do

! Close the file

500 close(iunit)

end subroutine airfoil_points

!=============================================================================80
!
! Subroutine to read an airfoil.  Assumes the number of points is already known.
! Also checks for incorrect format.
!
!=============================================================================80
subroutine airfoil_read(filename, npoints, labeled, x, z)

  character(*), intent(in) :: filename
  integer, intent(in) :: npoints
  logical, intent(in) :: labeled
  double precision, dimension(:), intent(inout) :: x, z

  integer :: i, iunit, ioerr, nswitch
  double precision :: dir1, dir2

! Open airfoil file

  iunit = 12
  open(unit=iunit, file=filename, status='old', position='rewind', iostat=ioerr)
  if (ioerr /= 0) then
     write(*,*) 'Error: cannot find airfoil file '//trim(filename)
     write(*,*)
     stop
  end if

! Read points from file

  if (labeled) read(iunit,*)
  do i = 1, npoints
    read(iunit,*,end=500,err=500) x(i), z(i)
  end do

! Close file

  close(iunit)

! Check that coordinates are formatted in a loop

  nswitch = 0
  dir1 = x(2) - x(1)
  do i = 3, npoints
    dir2 = x(i) - x(i-1)
    if (dir2 /= 0.d0) then
      if (dir2*dir1 < 0.d0) nswitch = nswitch + 1
      dir1 = dir2
    end if
  end do

  if (nswitch /= 1) then
!   Open the file again only to avoid error at label 500.
    open(unit=iunit, file=filename, status='old')
  else
    return
  end if

500 close(iunit)
  write(*,'(A)') "Error: incorrect format in "//trim(filename)//". File should"
  write(*,'(A)') "have x and y coordinates in 2 columns to form a single loop,"
  write(*,'(A)') "and there should be no blank lines.  See the user guide for"
  write(*,'(A)') "more information."
  stop

end subroutine airfoil_read

!=============================================================================80
!
! Changes airfoil point ordering to counterclockwise if necessary
!
!=============================================================================80
subroutine cc_ordering(foil)

  use vardef, only : airfoil_type

  type(airfoil_type), intent(inout) :: foil

  double precision, dimension(foil%npoint) :: xtemp, ztemp
  integer :: i, npoints

  npoints = foil%npoint

! Check if ordering needs to be switched

  if (foil%z(2) < foil%z(npoints-1)) then
    
    write(*,*) 'Changing point ordering to counter-clockwise ...'
    
    xtemp = foil%x
    ztemp = foil%z
    do i = 1, npoints
      foil%x(i) = xtemp(npoints-i+1)
      foil%z(i) = ztemp(npoints-i+1)
    end do

  end if

end subroutine cc_ordering

!=============================================================================80
!
! Subroutine to find leading edge of airfoil
!
! Input: airfoil(X,Z)
! Output: le: index of point closest to leading edge
!         xle: x-location of leading edge
!         zle: z-location of leading edge
!         addpoint_loc: integer giving the position at which to add a new point
!            for the leading edge. +1 means the index after le, -1 means the
!            index before le, and 0 means no new point is needed (x(le), z(le) 
!            is exactly at the leading edge).
!
! Note: addpoint_loc calculation assumes that airfoil points are ordered
! counter-clockwise
!
!=============================================================================80
subroutine le_find(x, z, le, xle, zle, addpoint_loc)

  use math_deps,    only : lmult

  double precision, dimension(:), intent(in) :: x, z
  integer, intent(out) :: le, addpoint_loc
  double precision, intent(out) :: xle, zle

  integer i, j, zeroval
  logical switch
  double precision, dimension(6) :: xp, zp
  double precision, dimension(6,6) :: A
  double precision, dimension(6) :: C
  double precision, dimension(2) :: bounds

  switch = .false.

! Find point where x-location of airfoil point starts increasing

  i = 1
  do while (.not. switch)
     if (x(i+1) >= x(i)) then

        le = i
        switch = .true.

!       Fit 5th-order polynomial to leading edge

        xp(1) = x(le+3)
        zp(1) = z(le+3)
        xp(2) = x(le+2)
        zp(2) = z(le+2)
        xp(3) = x(le+1)
        zp(3) = z(le+1)
        xp(4) = x(le)
        zp(4) = z(le)
        xp(5) = x(le-1)
        zp(5) = z(le-1)
        xp(6) = x(le-2)
        zp(6) = z(le-2)

!       Check for zero values - this row will get moved to bottom

        zeroval = 7
        do j = 1, 6
          if (zp(j) == 0.d0) then
            zeroval = j
            exit
          end if
        end do

!       Set up the system of equations. If there is a row with z(j) = 0, it
!       gets moved to the bottom.

        do j = 1, 6
          if (j - zeroval < 0) then
            A(j,:) = (/ zp(j)**5.d0, zp(j)**4.d0, zp(j)**3.d0, zp(j)**2.d0,    &
                        zp(j), 1.d0 /)
            C(j) = xp(j)
          elseif (j - zeroval > 0) then
            A(j-1,:) = (/ zp(j)**5.d0, zp(j)**4.d0, zp(j)**3.d0, zp(j)**2.d0,  &
                          zp(j), 1.d0 /)
            C(j-1) = xp(j)
          else
            A(6,:) = (/ zp(j)**5.d0, zp(j)**4.d0, zp(j)**3.d0, zp(j)**2.d0,    &
                        zp(j), 1.d0 /)
            C(6) = xp(j)
          end if
        end do

        polynomial_coefs = lmult(A,C)

!       Find minimum of polynomial using golden search

        bounds = (/ zp(1), zp(6) /)
        call golden_search(quintic, bounds, zle, xle)

!       Determine whether a new point needs to be added for the leading edge

        if (zle > z(le)) then
          addpoint_loc = -1
        elseif (zle == z(le)) then
          addpoint_loc = 0
        else
          addpoint_loc = 1
        end if

     else
        i = i + 1
     end if
  end do

end subroutine le_find

!=============================================================================80
!
! Golden search: bounded optimization in one dimension
!
! Inputs: 
!   objfunc: objective function
!   bounds: lower and upper bounds for the search
!
! Output: xmin, fmin
!
!=============================================================================80
subroutine golden_search(objfunc, bounds, xopt, fmin)

  double precision, dimension(2), intent(in) :: bounds
  double precision, intent(out) :: xopt, fmin

  interface
    double precision function objfunc(x)
      double precision, intent(in) :: x
    end function
  end interface

  double precision :: tol, T, dist
  double precision, dimension(4) :: xval, fval
  integer i, imax

  tol = 1.0e-09

  imax = 100
  dist = 1000.0d0

! Initialize search

  xval = (/ bounds(1), 0.d0, 0.d0, bounds(2) /)
  dist = bounds(2) - bounds(1)
  T = (3.d0 - sqrt(5.d0))/2.d0                  !Golden section
  xval(2) = (1.d0-T)*xval(1) + T*xval(4)
  xval(3) = T*xval(1) + (1.d0-T)*xval(4)
  fval = (/ objfunc(xval(1)), objfunc(xval(2)), objfunc(xval(3)),              &
            objfunc(xval(4)) /)
  if (fval(1) > fval(4)) then
     xopt = xval(4)
     fmin = fval(4)
  else
     xopt = xval(1)
     fmin = fval(1)
  end if

  i = 2
  do while (i < imax .and. dist > tol)

!    Eliminate the appropriate region

     if (fval(2) > fval(3)) then
        xopt = xval(3)
        fmin = fval(3)
        xval(1) = xval(2)
        xval(2) = xval(3)
        fval(1) = fval(2)
        fval(2) = fval(3)
        xval(3) = T*xval(1) + (1.d0-T)*xval(4)
        fval(3) = objfunc(xval(3))
     else
        xopt = xval(2)
        fmin = fval(2)
        xval(4) = xval(3)
        xval(3) = xval(2)
        fval(4) = fval(3)
        fval(3) = fval(2)
        xval(2) = (1.d0-T)*xval(1) + T*xval(4)
        fval(2) = objfunc(xval(2))
     end if
     dist = xval(4) - xval(1)

!    Print out warning if maximum number of iterations is reached

     if (i == imax .and. dist > tol) then
        write(*,*)
        write(*,*) 'Warning: golden search reached maximum number of iterations'
        write(*,*) '  without converging to the specified absolute error.'
        write(*,*)
     end if

     i = i + 1

  end do

end subroutine golden_search

!=============================================================================80
!
! Fifth order polynomial
!
!=============================================================================80
function quintic(x)

  double precision, intent(in) :: x
  double precision :: quintic

  quintic = polynomial_coefs(1)*x**5.d0 + polynomial_coefs(2)*x**4.d0 +        &
            polynomial_coefs(3)*x**3.d0 + polynomial_coefs(4)*x**2.d0 +        &
            polynomial_coefs(5)*x + polynomial_coefs(6)

end function quintic

!=============================================================================80
!
! Translates and scales an airfoil such that it has a length of 1 and the 
! leading edge is at the origin
!
!=============================================================================80
subroutine transform_airfoil(foil)

  use vardef, only : airfoil_type

  type(airfoil_type), intent(inout) :: foil

  integer :: npoints, i
  double precision :: scale_factor

  npoints = foil%npoint

! Translate so that the leading edge is at the origin

  do i = 1, npoints
    foil%x(i) = foil%x(i) - foil%xle
    foil%z(i) = foil%z(i) - foil%zle
  end do
  foil%xle = 0.d0
  foil%zle = 0.d0

! Scale airfoil so that it has a length of 1

  scale_factor = 1.d0 / maxval(foil%x)
  do i = 1, npoints
    foil%x(i) = foil%x(i)*scale_factor
    foil%z(i) = foil%z(i)*scale_factor
  end do

end subroutine transform_airfoil

!=============================================================================80
!
! Subroutine to determine the number of points on the top and bottom surfaces of
! an airfoil
!
!=============================================================================80
subroutine get_split_points(foil, pointst, pointsb, symmetrical)

  use vardef, only : airfoil_type

  type(airfoil_type), intent(in) :: foil
  integer, intent(out) :: pointst, pointsb
  logical, intent(in) :: symmetrical

  if (foil%addpoint_loc == 0) then
    pointst = foil%leclose
    pointsb = foil%npoint - foil%leclose + 1
  elseif (foil%addpoint_loc == -1) then
    pointst = foil%leclose 
    pointsb = foil%npoint - foil%leclose + 2
  else
    pointst = foil%leclose + 1
    pointsb = foil%npoint - foil%leclose + 1
  end if

! Modify for symmetrical airfoil (top surface will be mirrored)

  if (symmetrical) pointsb = pointst

end subroutine get_split_points

!=============================================================================80
!
! Subroutine to split an airfoil into top and bottom surfaces
!
!=============================================================================80
subroutine split_airfoil(foil, xseedt, xseedb, zseedt, zseedb, symmetrical)

  use vardef, only : airfoil_type

  type(airfoil_type), intent(in) :: foil
  double precision, dimension(:), intent(inout) :: xseedt, xseedb, zseedt,     &
                                                   zseedb
  logical, intent(in) :: symmetrical
  
  integer i, boundst, boundsb, pointst, pointsb

  pointst = size(xseedt,1)
  pointsb = size(xseedb,1)

  if (foil%addpoint_loc == 0) then
    boundst = foil%leclose - 1
    boundsb = foil%leclose + 1
  elseif (foil%addpoint_loc == -1) then
    boundst = foil%leclose - 1
    boundsb = foil%leclose
  else
    boundst = foil%leclose
    boundsb = foil%leclose + 1
  end if

! Copy points for the top surface

  xseedt(1) = foil%xle
  zseedt(1) = foil%zle
  do i = 1, pointst - 1
    xseedt(i+1) = foil%x(boundst-i+1)
    zseedt(i+1) = foil%z(boundst-i+1)
  end do

! Copy points for the bottom surface

  xseedb(1) = foil%xle
  zseedb(1) = foil%zle
  if (.not. symmetrical) then
    do i = 1, pointsb - 1
      xseedb(i+1) = foil%x(boundsb+i-1)
      zseedb(i+1) = foil%z(boundsb+i-1)
    end do
  else
    do i = 1, pointsb - 1
      xseedb(i+1) = xseedt(i+1)
      zseedb(i+1) = -zseedt(i+1)
    end do
  end if

end subroutine split_airfoil

!=============================================================================80
!
! Writes an airfoil to a labeled file
!
!=============================================================================80
subroutine airfoil_write(filename, title, foil)

  use vardef, only : airfoil_type

  character(*), intent(in) :: filename, title
  type(airfoil_type), intent(in) :: foil
  
  integer :: i, iunit

! Write notification to screen

  write(*,*)
  write(*,*) 'Writing labeled airfoil file '//trim(filename)//' ...'
  write(*,*)

! Open file for writing

  iunit = 13
  open(unit=iunit, file=filename, status='replace')

! Write label to file
  
  write(iunit,'(A)') trim(title)

! Write coordinates

  do i = 1, foil%npoint
    write(iunit,*) foil%x(i), foil%z(i)
  end do

! Close file

  close(iunit)

end subroutine airfoil_write

!=============================================================================80
!
! Allocates memory for buffer airfoil
!
!=============================================================================80
subroutine allocate_airfoil(foil)

  use vardef, only : airfoil_type

  type(airfoil_type), intent(inout) :: foil

  integer :: npoint
 
  npoint = foil%npoint
  allocate(foil%x(npoint))
  allocate(foil%z(npoint))

end subroutine allocate_airfoil

!=============================================================================80
!
! Deallocates memory for buffer airfoil
!
!=============================================================================80
subroutine deallocate_airfoil(foil)

  use vardef, only : airfoil_type

  type(airfoil_type), intent(inout) :: foil

  deallocate(foil%x)
  deallocate(foil%z)

end subroutine deallocate_airfoil

!=============================================================================80
!
! Checks if a given character is a number
!
!=============================================================================80
function isnum(s)

  character, intent(in) :: s
  logical :: isnum

  select case (s)
    case ('0')
      isnum = .true.
    case ('1')
      isnum = .true.
    case ('2')
      isnum = .true.
    case ('3')
      isnum = .true.
    case ('4')
      isnum = .true.
    case ('5')
      isnum = .true.
    case ('6')
      isnum = .true.
    case ('7')
      isnum = .true.
    case ('8')
      isnum = .true.
    case ('9')
      isnum = .true.
    case default
      isnum = .false.
  end select

end function isnum

!=============================================================================80
!
! Stops and prints an error message, or just warns
!
!=============================================================================80
subroutine my_stop(message, stoptype)

  character(*), intent(in) :: message
  character(4), intent(in), optional :: stoptype

  if ((.not. present(stoptype)) .or. (stoptype == 'stop')) then
    write(*,'(A)') 'Error: '//trim(message)
    write(*,*)
    stop
  else
    write(*,'(A)') 'Warning: '//trim(message)
    write(*,*)
  end if

end subroutine my_stop


end module airfoil_operations
