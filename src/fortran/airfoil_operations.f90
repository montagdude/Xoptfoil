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
subroutine get_seed_airfoil(seed_airfoil, airfoil_file, naca_options, foil,    &
                            xoffset, zoffset, foilscale)

  use vardef,       only : airfoil_type
  use xfoil_driver, only : smooth_paneling
  use naca,         only : naca_options_type, naca_456

  character(*), intent(in) :: seed_airfoil, airfoil_file
  type(naca_options_type), intent(in) :: naca_options
  type(airfoil_type), intent(out) :: foil
  double precision, intent(out) :: xoffset, zoffset, foilscale

  type(airfoil_type) :: tempfoil
  integer :: pointsmcl

  if (trim(seed_airfoil) == 'from_file') then

!   Read seed airfoil from file

    call load_airfoil(airfoil_file, tempfoil)

  elseif (trim(seed_airfoil) == 'naca') then

!   Create NACA 4, 4M, 5, 6, or 6A series airfoil

    pointsmcl = 200
    call naca_456(naca_options, pointsmcl, tempfoil)

  else

    write(*,*) "Error: seed_airfoil should be 'from_file' or 'naca'."
    write(*,*)
    stop

  end if

! Use Xfoil to smooth airfoil paneling

  call smooth_paneling(tempfoil, 200, foil)

! Calculate leading edge information

  call le_find(foil%x, foil%z, foil%leclose, foil%xle, foil%zle,               &
               foil%addpoint_loc)

! Translate and scale

  call transform_airfoil(foil, xoffset, zoffset, foilscale)

end subroutine get_seed_airfoil

!=============================================================================80
!
! Reads an airfoil from a file, loads it into the airfoil_type, sets ordering
! correctly
!
!=============================================================================80
subroutine load_airfoil(filename, foil)

  use vardef,      only : airfoil_type
  use memory_util, only : allocate_airfoil

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

  use vardef,    only : airfoil_type
  use math_deps, only : norm_2

  type(airfoil_type), intent(inout) :: foil

  double precision, dimension(foil%npoint) :: xtemp, ztemp
  double precision, dimension(2) :: tevec1, tevec2
  double precision :: len1, len2
  integer :: i, npoints

  npoints = foil%npoint

! Check if ordering needs to be switched

  tevec1(1) = foil%x(2) - foil%x(1)
  tevec1(2) = foil%z(2) - foil%z(1)
  len1 = norm_2(tevec1)

  tevec2(1) = foil%x(npoints-1) - foil%x(npoints)
  tevec2(2) = foil%z(npoints-1) - foil%z(npoints)
  len2 = norm_2(tevec2)

  if ( (len1 == 0.d0) .or. (len2 == 0.d0) )                                    &
    call my_stop("Panel with 0 length detected near trailing edge.")

  tevec1 = tevec1/len1
  tevec2 = tevec2/len2

  if (tevec1(2) < tevec2(2)) then
    
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
!=============================================================================80
subroutine le_find(x, z, le, xle, zle, addpoint_loc)

  use math_deps, only : norm_2

  double precision, dimension(:), intent(in) :: x, z
  integer, intent(out) :: le, addpoint_loc
  double precision, intent(out) :: xle, zle

  integer :: i, npt
  double precision, dimension(:), allocatable :: s, xp, zp
  double precision, dimension(2) :: r1, r2
  double precision :: sle, dist1, dist2, dot

! Get leading edge location from Xfoil

  npt = size(x,1)
  allocate(s(npt))
  allocate(xp(npt))
  allocate(zp(npt))
  call SCALC(x, z, s, npt)
  call SEGSPL(x, xp, s, npt)
  call SEGSPL(z, zp, s, npt)
  call LEFIND(sle, x, xp, z, zp, s, npt, .true.)
  call SEVAL_INTERFACE(xle, sle, x, xp, s, npt)
  call SEVAL_INTERFACE(zle, sle, z, zp, s, npt)
  deallocate(s)
  deallocate(xp)
  deallocate(zp)

! Determine leading edge index and where to add a point

  npt = size(x,1)
  do i = 1, npt-1
    r1(1) = xle - x(i)
    r1(2) = zle - z(i)
    dist1 = norm_2(r1)
    if (dist1 /= 0.d0) r1 = r1/dist1

    r2(1) = xle - x(i+1)
    r2(2) = zle - z(i+1)
    dist2 = norm_2(r2)
    if (dist2 /= 0.d0) r2 = r2/dist2

    dot = dot_product(r1, r2)

    if (dist1 == 0.d0) then
      le = i
      addpoint_loc = 0
      exit
    else if (dist2 == 0.d0) then
      le = i+1
      addpoint_loc = 0
      exit
    else if (dot < 0.d0) then
      if (dist1 < dist2) then
        le = i
        addpoint_loc = 1
      else
        le = i+1
        addpoint_loc = -1
      end if
      exit
    end if
  end do

end subroutine le_find

!=============================================================================80
!
! Translates and scales an airfoil such that it has a length of 1 and the 
! leading edge is at the origin. Also outputs transformations performed.
!
!=============================================================================80
subroutine transform_airfoil(foil, xoffset, zoffset, foilscale)

  use vardef, only : airfoil_type

  type(airfoil_type), intent(inout) :: foil
  double precision, intent(out) :: xoffset, zoffset, foilscale

  integer :: npoints, i

  npoints = foil%npoint

! Translate so that the leading edge is at the origin

  do i = 1, npoints
    foil%x(i) = foil%x(i) - foil%xle
    foil%z(i) = foil%z(i) - foil%zle
  end do
  xoffset = -foil%xle
  zoffset = -foil%zle
  foil%xle = 0.d0
  foil%zle = 0.d0

! Scale airfoil so that it has a length of 1

  foilscale = 1.d0 / maxval(foil%x)
  do i = 1, npoints
    foil%x(i) = foil%x(i)*foilscale
    foil%z(i) = foil%z(i)*foilscale
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
    write(*,*)
    write(*,'(A)') 'Error: '//trim(message)
    write(*,*)
    stop
  else
    write(*,*)
    write(*,'(A)') 'Warning: '//trim(message)
    write(*,*)
  end if

end subroutine my_stop

end module airfoil_operations
