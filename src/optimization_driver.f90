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

module optimization_driver

! Contains subroutines to set options and conditions for optimization and to
! issue optimizer calls

  implicit none

  contains

!=============================================================================80
!
! Subroutine to drive the optimization
!
!=============================================================================80
subroutine optimize(search_type, global_search, local_search, matchfoil_file,  &
                    initial_perturb, pso_options, ds_options, optdesign, fmin, &
                    steps, fevals)

  use vardef,             only : airfoil_type, match_foils, xmatcht, xmatchb,  &
                                 zmatcht, zmatchb, xseedt, xseedb, zseedt,     &
                                 zseedb, growth_allowed, shape_functions,      &
                                 symmetrical
  use optimization,       only : pso_options_type, ds_options_type,            &
                                 particleswarm, simplex_search
  use airfoil_evaluation, only : objective_function
  use airfoil_operations, only : get_seed_airfoil, get_split_points,           &
                                 split_airfoil, allocate_airfoil,              &
                                 deallocate_airfoil, my_stop
  use math_deps,          only : interp_vector
  use input_sanity,       only : check_seed

  character(*), intent(in) :: search_type, global_search, local_search,        &
                              matchfoil_file
  double precision, intent(in) :: initial_perturb
  type(pso_options_type), intent(in) :: pso_options
  type(ds_options_type), intent(in) :: ds_options
  double precision, dimension(:), intent(inout) :: optdesign
  double precision, intent(out) :: fmin
  integer, intent(out) :: steps, fevals

  type(airfoil_type) :: match_foil
  integer :: pointst, pointsb, counter, nfuncs
  double precision, dimension(:), allocatable :: zttmp, zbtmp
  double precision, dimension(size(optdesign,1)) :: xmin, xmax, x0
  double precision :: len1, len2, growth1, growth2, t1fact, t2fact
  integer :: stepsg, fevalsg, stepsl, fevalsl, i

!  Preliminary things for non-aerodynamic optimization

  if (match_foils) then

    write(*,*) 'Note: using the optimizer to match the seed airfoil to the '
    write(*,*) 'airfoil about to be loaded.'
    write(*,*)

!   Check if symmetrical airfoil was requested (not allowed for this type)

    if (symmetrical)                                                           &
      call my_stop("Symmetrical airfoil constraint not permitted for non-"//&
                   "aerodynamic optimizations.")

!   Load airfoil to match

    call get_seed_airfoil('from_file', matchfoil_file, '0000', match_foil)

!   Split match_foil into upper and lower halves

    call get_split_points(match_foil, pointst, pointsb, .false.)
    allocate(xmatcht(pointst))
    allocate(zmatcht(pointst))
    allocate(xmatchb(pointsb))
    allocate(zmatchb(pointsb))
    call split_airfoil(match_foil, xmatcht, xmatchb, zmatcht, zmatchb, .false.)

!   Deallocate derived type version of airfoil being matched

    call deallocate_airfoil(match_foil)

!   Interpolate x-vals of foil to match to seed airfoil points to x-vals

    pointst = size(xseedt,1)
    pointsb = size(xseedb,1)
    allocate(zttmp(pointst))
    allocate(zbtmp(pointsb))
    zttmp(pointst) = zmatcht(size(zmatcht,1))
    zbtmp(pointsb) = zmatchb(size(zmatchb,1))
    call interp_vector(xmatcht, zmatcht, xseedt(1:pointst-1),                  &
                       zttmp(1:pointst-1))
    call interp_vector(xmatchb, zmatchb, xseedb(1:pointsb-1),                  &
                       zbtmp(1:pointsb-1))

!   Re-set coordinates of foil to match from interpolated points
    
    deallocate(xmatcht)
    deallocate(zmatcht)
    deallocate(xmatchb)
    deallocate(zmatchb)
    allocate(xmatcht(pointst))
    allocate(zmatcht(pointst))
    allocate(xmatchb(pointsb))
    allocate(zmatchb(pointsb))
    xmatcht = xseedt
    xmatchb = xseedb
    zmatcht = zttmp
    zmatchb = zbtmp

!   Deallocate temporary arrays

    deallocate(zttmp)
    deallocate(zbtmp)

  else

    write(*,*) "Optimizing for requested operating points."
    write(*,*)

!   Get allowable panel growth rate

    pointst = size(xseedt,1)
    pointsb = size(xseedb,1)
    growth_allowed = 0.d0

!   Top surface growth rates

    len1 = sqrt((xseedt(2)-xseedt(1))**2.d0 + (zseedt(2)-zseedt(1))**2.d0)
    do i = 2, pointst - 1
      len2 = sqrt((xseedt(i+1)-xseedt(i))**2.d0 + (zseedt(i+1)-zseedt(i))**2.d0)
      growth1 = len2/len1
      growth2 = len1/len2
      if (max(growth1,growth2) > growth_allowed)                               &
          growth_allowed = 1.5d0*max(growth1,growth2)
      len1 = len2
    end do

!   Bottom surface growth rates

    len1 = sqrt((xseedb(2)-xseedb(1))**2.d0 + (zseedb(2)-zseedb(1))**2.d0)
    do i = 2, pointsb - 1
      len2 = sqrt((xseedb(i+1)-xseedb(i))**2.d0 + (zseedb(i+1)-zseedb(i))**2.d0)
      growth1 = len2/len1
      growth2 = len1/len2
      if (max(growth1,growth2) > growth_allowed)                               &
          growth_allowed = 1.5d0*max(growth1,growth2)
      len1 = len2
    end do

  end if

! Make sure seed airfoil passes constraints, and get scaling factors for
! operating points

  call check_seed()

! Perform optimization: global, local, or global + local

  stepsg = 0
  fevalsg = 0
  stepsl = 0
  fevalsl = 0

  if (trim(search_type) == 'global_and_local' .or. trim(search_type) ==        &
      'global') then

    if (trim(global_search) == 'particle_swarm') then

!     Set up mins and maxes
      
      if (trim(shape_functions) == 'naca') then

        xmin(:) = -initial_perturb/2.d0
        xmax(:) = initial_perturb/2.d0

      else

        nfuncs = size(optdesign,1)/3

!       Scaling factors for Hicks-Henne bump locations and widths

        t1fact = 0.05d0
        t2fact = 0.005d0

        do i = 1, nfuncs
          counter = 3*(i-1)
          xmin(counter+1) = -initial_perturb/2.d0
          xmax(counter+1) = initial_perturb/2.d0
          xmin(counter+2) = 0.0001d0*t1fact
          xmax(counter+2) = 1.d0*t1fact
          xmin(counter+3) = 0.05d0*t2fact
          xmax(counter+3) = 10.d0*t2fact
        end do

      end if

!     Particle swarm optimization

      call particleswarm(optdesign, fmin, stepsg, fevalsg, objective_function, &
                         xmin, xmax, pso_options)

    end if

  end if

  if (trim(search_type) == 'global_and_local' .or. trim(search_type) ==        &
      'local') then

    if (trim(local_search) == 'simplex') then

!     Simplex optimization

      if (trim(search_type) == 'global_and_local') then
        x0 = optdesign
      else

        nfuncs = size(optdesign,1)/3

!       Scaling factors for Hicks-Henne bump locations and widths

        t1fact = 0.05d0
        t2fact = 0.005d0

        do i = 1, nfuncs
          counter = 3*(i-1)
          x0(counter+1) = 0.d0
          x0(counter+2) = 0.5d0*t1fact
          x0(counter+2) = 1.d0*t2fact
        end do

      end if

      call simplex_search(optdesign, fmin, stepsl, fevalsl, objective_function,&
                          x0, ds_options)

    end if

  end if

! Total number of steps and function evaluations

  steps = stepsg + stepsl
  fevals = fevalsg + fevalsl

! Deallocate memory for match_foil

  if (match_foils) then
    deallocate(xmatcht)
    deallocate(zmatcht)
    deallocate(xmatchb)
    deallocate(zmatchb)
  end if

end subroutine optimize

!=============================================================================80
!
! Writes final airfoil design to a file 
!
!=============================================================================80
subroutine write_final_design(optdesign, shapetype, output_prefix, symmetrical)

  use vardef,             only : airfoil_type, xseedt, xseedb, zseedt, zseedb
  use airfoil_operations, only : allocate_airfoil, deallocate_airfoil,         &
                                 airfoil_write
  use parameterization,   only : top_shape_function, bot_shape_function,       &
                                 create_airfoil

  double precision, dimension(:), intent(in) :: optdesign
  character(*), intent(in) :: shapetype, output_prefix
  logical, intent(in) :: symmetrical

  double precision, dimension(size(xseedt,1)) :: zt_new
  double precision, dimension(size(xseedb,1)) :: zb_new
  integer :: dvtbnd1, dvtbnd2, dvbbnd1, dvbbnd2, nmodest, nmodesb, nptt, nptb, i
  type(airfoil_type) :: final_airfoil
  character(80) :: output_file

  nmodest = size(top_shape_function,1)
  nmodesb = size(bot_shape_function,1)
  nptt = size(xseedt,1)
  nptb = size(xseedb,1)

! Set modes for top and bottom surfaces

  dvtbnd1 = 1
  if (trim(shapetype) == 'naca') then
    dvtbnd2 = nmodest
    dvbbnd2 = nmodest + nmodesb
  else
    dvtbnd2 = nmodest*3
    dvbbnd2 = nmodest*3 + nmodesb*3
  end if
  dvbbnd1 = dvtbnd2 + 1

! Overwrite lower DVs for symmetrical airfoils (they are not used)

  if (symmetrical) then
    dvbbnd1 = 1
    dvbbnd2 = dvtbnd2
  end if

! Create top and bottom surfaces by perturbation of seed airfoil

  call create_airfoil(xseedt, zseedt, xseedb, zseedb,                          &
                      optdesign(dvtbnd1:dvtbnd2), optdesign(dvbbnd1:dvbbnd2),  &
                      zt_new, zb_new, shapetype, symmetrical)

! Format coordinates in a single loop (in airfoil_type derived type)

  final_airfoil%npoint = nptt + nptb - 1
  call allocate_airfoil(final_airfoil)
  do i = 1, nptt
    final_airfoil%x(i) = xseedt(nptt-i+1)
    final_airfoil%z(i) = zt_new(nptt-i+1)
  end do
  do i = 1, nptb - 1
   final_airfoil%x(i+nptt) = xseedb(i+1)
   final_airfoil%z(i+nptt) = zb_new(i+1)
  end do

! Write to file

  output_file = trim(output_prefix)//'.dat'
  call airfoil_write(output_file, output_prefix, final_airfoil)

! Deallocate final airfoil

  call deallocate_airfoil(final_airfoil)

end subroutine write_final_design

end module optimization_driver
