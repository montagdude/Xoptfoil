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

module input_sanity

  implicit none

  contains

!=============================================================================80
!
! Checks that the seed airfoil passes all constraints, sets scale factors for
! objective functions at each operating point, and optionally sets the
! minimum allowable pitching moment.
!
!=============================================================================80
subroutine check_seed()

  use vardef
  use math_deps,          only : interp_vector, curvature
  use xfoil_driver,       only : xfoil_init, xfoil_cleanup, run_xfoil
  use xfoil_inc,          only : AMAX
  use airfoil_operations, only : allocate_airfoil, deallocate_airfoil, my_stop
  use airfoil_evaluation, only : xfoil_options, xfoil_geom_options

  double precision, dimension(:), allocatable :: x_interp, thickness
  double precision, dimension(:), allocatable :: zt_interp, zb_interp
  double precision, dimension(size(xseedt,1)+size(xseedb,1)-1) :: curv
  double precision :: penaltyval, tegap, gapallow, maxthick, heightfactor
  double precision :: panang1, panang2, maxpanang, curv1, curv2
  double precision :: checkval
  double precision, dimension(noppoint) :: lift, drag, moment, viscrms
  integer :: i, nptt, nptb, nreversals, nptint
  character(30) :: text
  character(4) :: stoptype

  stoptype = seed_violation_handling
  penaltyval = 0.d0
  nptt = size(xseedt,1)
  nptb = size(xseedb,1)

  write(*,*) 'Checking to make sure seed airfoil passes all constraints ...'
  write(*,*)

! Allocate memory needed for check

  curr_foil%npoint = size(xseedt,1) + size(xseedb,1) - 1
  call allocate_airfoil(curr_foil)
  call xfoil_init()

! Format coordinates in a single loop in derived type

  do i = 1, nptt
    curr_foil%x(i) = xseedt(nptt-i+1)
    curr_foil%z(i) = zseedt(nptt-i+1)
  end do
  do i = 1, nptb-1
    curr_foil%x(i+nptt) = xseedb(i+1)
    curr_foil%z(i+nptt) = zseedb(i+1)
  end do

! Too blunt or sharp leading edge

  panang1 = atan((zseedt(2)-zseedt(1))/(xseedt(2)-xseedt(1))) *                &
            180.d0/acos(-1.d0)
  panang2 = atan((zseedb(1)-zseedb(2))/(xseedb(2)-xseedb(1))) *                &
            180.d0/acos(-1.d0)
  maxpanang = max(panang2,panang1)
  if (maxpanang > 89.99d0)                                                     &
    call my_stop("Seed airfoil's leading edge is too blunt.", stoptype)
  if (abs(panang1 - panang2) > 20.d0)                                          &
    call my_stop("Seed airfoil's leading edge is too sharp.", stoptype)

! Interpolate either bottom surface to top surface x locations or vice versa
! to determine thickness

  if (xseedt(nptt) <= xseedb(nptb)) then
    allocate(x_interp(nptt))
    allocate(zt_interp(nptt))
    allocate(zb_interp(nptt))
    allocate(thickness(nptt))
    nptint = nptt
    call interp_vector(xseedb, zseedb, xseedt, zb_interp)
    x_interp = xseedt
    zt_interp = zseedt
  else
    allocate(x_interp(nptb))
    allocate(zt_interp(nptb))
    allocate(zb_interp(nptb))
    allocate(thickness(nptb))
    nptint = nptb
    call interp_vector(xseedt, zseedt, xseedb, zt_interp)
    x_interp = xseedb
    zb_interp = zseedb
  end if

! Compute thickness parameters

  tegap = zseedt(nptt) - zseedb(nptb)
  maxthick = 0.d0
  heightfactor = tan(min_te_angle*acos(-1.d0)/180.d0/2.d0)

  do i = 2, nptint - 1

!   Thickness array and max thickness
    
    thickness(i) = zt_interp(i) - zb_interp(i)
    if (thickness(i) > maxthick) maxthick = thickness(i)

!   Check if thinner than specified wedge angle on back half of airfoil
    
    if (x_interp(i) > 0.5d0) then
      gapallow = tegap + 2.d0 * heightfactor * (x_interp(nptint) -             &
                                                x_interp(i))
      if (thickness(i) < gapallow)                                             &
        penaltyval = penaltyval + (gapallow - thickness(i))/0.001d0
    end if

  end do

! Free memory

  deallocate(x_interp)
  deallocate(zt_interp)
  deallocate(zb_interp)
  deallocate(thickness)

! Too thin on back half

  if (penaltyval > 0.d0)                                                       &
     call my_stop("Seed airfoil is thinner than min_te_angle near the "//      &
                  "trailing edge.", stoptype)
  penaltyval = 0.d0

! Max thickness too low

  if (maxthick < min_thickness)                                                &
    call my_stop("Seed airfoil violates min_thickness constraint.", stoptype)

! Max thickness too high

  if (maxthick > max_thickness)                                                &
    call my_stop("Seed airfoil violates max_thickness constraint.", stoptype)

! Check for curvature reversals

  if (check_curvature) then

!   Compute curvature

    curv = curvature(curr_foil%npoint, curr_foil%x, curr_foil%z)

!   Check number of reversals that exceed the threshold

    nreversals = 0
    curv1 = 0.d0

    do i = 2, nptt + nptb - 2

      if (abs(curv(i)) >= curv_threshold) then
        curv2 = curv(i)
        if (curv2*curv1 < 0.d0) nreversals = nreversals + 1
        curv1 = curv2
      end if

    end do

    if (nreversals > max_curv_reverse)                                         &
      call my_stop("Seed airfoil violates max_curv_reverse constraint.",       &
                   stoptype)

  end if

! Check for bad combinations of operating conditions and optimization types

  do i = 1, noppoint
    write(text,*) i
    text = adjustl(text)

    if ((op_point(i) == 0.d0) .and. (op_mode(i) == 'spec-cl') .and.            &
        (trim(optimization_type(i)) /= 'min-drag')) then
      write(*,*) "Error: operating points "//trim(text)//" is at Cl = 0. "//   &
                 "Need to use 'min-drag' optimization type in this case."
      write(*,*) 
      stop
    elseif ((op_mode(i) == 'spec-cl') .and.                                    &
            (trim(optimization_type(i)) == 'max-lift')) then
      write(*,*) "Error: Cl is specified for operating point "//trim(text)//   &
                 ". Cannot use 'max-lift' optimization type in this case."
      write(*,*) 
      stop
    end if

  end do

! Analyze airfoil at requested operating conditions with Xfoil

  call run_xfoil(curr_foil, xfoil_geom_options, op_point(1:noppoint),          &
                 op_mode(1:noppoint), reynolds(1:noppoint), mach(1:noppoint),  &
                 xfoil_options, lift, drag, moment, viscrms)

! Check for unconverged points

  do i = 1, noppoint
    if (viscrms(i) > 1.0D-04) then
      write(text,*) i
      text = adjustl(text)
      call my_stop("Xfoil calculations did not converge for operating point "//&
                   trim(text)//".", stoptype)
    end if
  end do

! Set moment constraint or check for violation of specified constraint

  if (trim(moment_constraint_type) == 'use_seed') then
    min_moment = minval(moment)
  elseif (trim(moment_constraint_type) == 'specify') then
    if (minval(moment) < min_moment)                                           &
      call my_stop("Seed airfoil violates min_moment constraint.", stoptype)
  end if

! Evaluate objectives to establish scale factors for each point

  do i = 1, noppoint
    if (trim(optimization_type(i)) == 'min-sink') then
      if (lift(i) > 0.d0) then
        checkval = drag(i)/lift(i)**1.5d0
      else
        write(text,*) i
        text = adjustl(text)
        write(*,*) "Error: operating point "//trim(text)//" has Cl <= 0. "//   &
                   "Cannot use min-sink optimization in this case."
        write(*,*)
        stop
      end if
    elseif (trim(optimization_type(i)) == 'max-glide') then
      checkval = drag(i)/lift(i)
    elseif (trim(optimization_type(i)) == 'min-drag') then
      checkval = drag(i)
    elseif (trim(optimization_type(i)) == 'max-lift') then
      checkval = 1.d0/lift(i)
    else
      write(*,*)
      write(*,*) "Error: requested optimization_type for operating point "//   &
                 trim(text)//" not recognized."
      stop
    end if
    scale_factor(i) = 1.d0/checkval
  end do

! Penalty for too large panel angles

  maxpanang = AMAX
  if (maxpanang > 25.d0)                                                       &
    call my_stop("Seed airfoil panel angles are too large. Try adjusting "//   &
                 "xfoil_paneling_options.", stoptype)

! Deallocate memory

  call deallocate_airfoil(curr_foil)
  call xfoil_cleanup()
  
end subroutine check_seed

end module input_sanity
