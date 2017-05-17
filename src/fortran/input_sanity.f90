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
subroutine check_seed(xoffset, zoffset, foilscale)

  use vardef
  use math_deps,          only : interp_vector, curvature, derv1f1, derv1b1
  use xfoil_driver,       only : run_xfoil
  use xfoil_inc,          only : AMAX, CAMBR
  use airfoil_evaluation, only : xfoil_options, xfoil_geom_options

  double precision, intent(in) :: xoffset, zoffset, foilscale

  double precision, dimension(:), allocatable :: x_interp, thickness
  double precision, dimension(:), allocatable :: zt_interp, zb_interp
  double precision, dimension(size(xseedt,1)) :: curvt
  double precision, dimension(size(xseedb,1)) :: curvb
  double precision, dimension(naddthickconst) :: add_thickvec
  double precision :: penaltyval, tegap, gapallow, maxthick, heightfactor
  double precision :: panang1, panang2, maxpanang, curv1, curv2
  double precision :: checkval, len1, len2, growth1, growth2, xtrans, ztrans
  double precision, dimension(noppoint) :: lift, drag, moment, viscrms, alpha, &
                                           xtrt, xtrb
  double precision :: pi
  integer :: i, nptt, nptb, nreversalst, nreversalsb, nptint
  character(30) :: text, text2
  character(14) :: opt_type
  logical :: addthick_violation

  penaltyval = 0.d0
  pi = acos(-1.d0)
  nptt = size(xseedt,1)
  nptb = size(xseedb,1)

  write(*,*) 'Checking to make sure seed airfoil passes all constraints ...'
  write(*,*)

! Get allowable panel growth rate

  growth_allowed = 0.d0

! Top surface growth rates

  len1 = sqrt((xseedt(2)-xseedt(1))**2.d0 + (zseedt(2)-zseedt(1))**2.d0)
  do i = 2, nptt - 1
    len2 = sqrt((xseedt(i+1)-xseedt(i))**2.d0 + (zseedt(i+1)-zseedt(i))**2.d0)
    growth1 = len2/len1
    growth2 = len1/len2
    if (max(growth1,growth2) > growth_allowed)                                 &
        growth_allowed = 1.5d0*max(growth1,growth2)
    len1 = len2
  end do

! Bottom surface growth rates

  len1 = sqrt((xseedb(2)-xseedb(1))**2.d0 + (zseedb(2)-zseedb(1))**2.d0)
  do i = 2, nptb - 1
    len2 = sqrt((xseedb(i+1)-xseedb(i))**2.d0 + (zseedb(i+1)-zseedb(i))**2.d0)
    growth1 = len2/len1
    growth2 = len1/len2
    if (max(growth1,growth2) > growth_allowed)                               &
        growth_allowed = 1.5d0*max(growth1,growth2)
    len1 = len2
  end do

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
  if (maxpanang > 89.99d0) then
    write(text,'(F8.4)') maxpanang
    text = adjustl(text)
    write(*,*) "LE panel angle: "//trim(text)//" degrees"
    call ask_stop("Seed airfoil's leading edge is too blunt.")
  end if
  if (abs(panang1 - panang2) > 20.d0) then
    write(text,'(F8.4)') abs(panang1 - panang2)
    text = adjustl(text)
    write(*,*) "LE panel angle: "//trim(text)//" degrees"
    call ask_stop("Seed airfoil's leading edge is too sharp.")
  end if

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
      if (thickness(i) < gapallow) then
        xtrans = x_interp(i)/foilscale - xoffset
        write(text,'(F8.4)') xtrans
        text = adjustl(text)
        write(*,*) "Detected too thin at x = "//trim(text)
        penaltyval = penaltyval + (gapallow - thickness(i))/0.1d0
      end if
    end if

  end do

! Too thin on back half

  if (penaltyval > 0.d0)                                                       &
     call ask_stop("Seed airfoil is thinner than min_te_angle near the "//&
                   "trailing edge.")
  penaltyval = 0.d0

! Check additional thickness constraints

  if (naddthickconst > 0) then
    call interp_vector(x_interp, thickness,                                    &
                       addthick_x(1:naddthickconst), add_thickvec)

    addthick_violation = .false.
    do i = 1, naddthickconst
      if ( (add_thickvec(i) < addthick_min(i)) .or.                            &
           (add_thickvec(i) > addthick_max(i)) ) then
        addthick_violation = .true.
        write(text,'(F8.4)') addthick_x(i)
        text = adjustl(text)
        write(text2,'(F8.4)') add_thickvec(i)
        text2 = adjustl(text2)
        write(*,*) "Thickness at x = "//trim(text)//": "//trim(text2)
      end if
    end do

    if (addthick_violation)                                                    &
      call ask_stop("Seed airfoil violates one or more thickness constraints.")
  end if
  
! Free memory

  deallocate(x_interp)
  deallocate(zt_interp)
  deallocate(zb_interp)
  deallocate(thickness)

! Max thickness too low

  if (maxthick < min_thickness) then
    write(text,'(F8.4)') maxthick
    text = adjustl(text)
    write(*,*) "Thickness: "//trim(text)
    call ask_stop("Seed airfoil violates min_thickness constraint.")
  end if

! Max thickness too high

  if (maxthick > max_thickness) then
    write(text,'(F8.4)') maxthick
    text = adjustl(text)
    write(*,*) "Thickness: "//trim(text)
    call ask_stop("Seed airfoil violates max_thickness constraint.")
  end if

! Check for curvature reversals

  if (check_curvature) then

!   Compute curvature on top and bottom surfaces

    curvt = curvature(nptt, xseedt, zseedt)
    curvb = curvature(nptb, xseedb, zseedb)

!   Check number of reversals that exceed the threshold

    nreversalst = 0
    curv1 = 0.d0
    do i = 2, nptt - 1
      if (abs(curvt(i)) >= curv_threshold) then
        curv2 = curvt(i)
        if (curv2*curv1 < 0.d0) then
          xtrans = xseedt(i)/foilscale - xoffset
          write(text,'(F8.4)') xtrans
          text = adjustl(text)
          ztrans = zseedt(i)/foilscale - zoffset
          write(text2,'(F8.4)') ztrans
          text2 = adjustl(text2)
          write(*,*) "Curvature reversal on top surface near (x, z) = ("//&
                         trim(text)//", "//trim(text2)//")"
          write(text,'(F8.4)') curvt(i)
          text = adjustl(text)
          write(*,*) "Curvature: "//trim(text)
          nreversalst = nreversalst + 1
        end if
        curv1 = curv2
      end if
    end do

    if (nreversalst > max_curv_reverse_top)                                    &
      call ask_stop("Seed airfoil violates max_curv_reverse_top constraint.")

!   Bottom surface

    nreversalsb = 0
    curv1 = 0.d0
    do i = 2, nptb - 1
      if (abs(curvb(i)) >= curv_threshold) then
        curv2 = curvb(i)
        if (curv2*curv1 < 0.d0) then
          xtrans = xseedb(i)/foilscale - xoffset
          write(text,'(F8.4)') xtrans
          text = adjustl(text)
          ztrans = zseedb(i)/foilscale - zoffset
          write(text2,'(F8.4)') ztrans
          text2 = adjustl(text2)
          write(*,*) "Curvature reversal on bot surface near (x, z) = ("//&
                         trim(text)//", "//trim(text2)//")"
          write(text,'(F8.4)') curvb(i)
          text = adjustl(text)
          write(*,*) "Curvature: "//trim(text)
          nreversalsb = nreversalsb + 1
        end if
        curv1 = curv2
      end if
    end do

    if (nreversalsb > max_curv_reverse_bot)                                    &
      call ask_stop("Seed airfoil violates max_curv_reverse_bot constraint.")

  end if

! Check for bad combinations of operating conditions and optimization types

  do i = 1, noppoint
    write(text,*) i
    text = adjustl(text)

    opt_type = optimization_type(i)
    if ((op_point(i) <= 0.d0) .and. (op_mode(i) == 'spec-cl')) then
      if ( (trim(opt_type) /= 'min-drag') .and.                                &
           (trim(opt_type) /= 'max-xtr') .and.                                 &
           (trim(opt_type) /= 'max-lift-slope') ) then
        write(*,*) "Error: operating point "//trim(text)//" is at Cl = 0. "//  &
                 "Cannot use '"//trim(opt_type)//"' optimization in this case."
        write(*,*) 
        stop
      end if
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
                 use_flap, x_flap, y_flap, y_flap_spec,                        &
                 flap_degrees(1:noppoint), xfoil_options, lift, drag, moment,  &
                 viscrms, alpha, xtrt, xtrb, ncrit_pt)

! Penalty for too large panel angles

  if (AMAX > 25.d0) then
    write(text,'(F8.4)') AMAX
    text = adjustl(text)
    write(*,*) "Max panel angle: "//trim(text)
    call ask_stop("Seed airfoil panel angles are too large. Try adjusting "//&
                  "xfoil_paneling_options.")
  end if

! Camber too high

  if (CAMBR > max_camber) then
    write(text,'(F8.4)') CAMBR
    text = adjustl(text)
    write(*,*) "Camber: "//trim(text)
    call ask_stop("Seed airfoil violates max_camber constraint.")
  end if

! Camber too low

  if (CAMBR < min_camber) then
    write(text,'(F8.4)') CAMBR
    text = adjustl(text)
    write(*,*) "Camber: "//trim(text)
    call ask_stop("Seed airfoil violates min_camber constraint.")
  end if

! Check for unconverged points

  do i = 1, noppoint
    if (viscrms(i) > 1.0D-04) then
      write(text,*) i
      text = adjustl(text)
      call ask_stop("Xfoil calculations did not converge for operating "//&
                    "point "//trim(text)//".")
    end if
  end do

! Set moment constraint or check for violation of specified constraint

  do i = 1, noppoint
    if (trim(moment_constraint_type(i)) == 'use_seed') then
      min_moment(i) = moment(i)
    elseif (trim(moment_constraint_type(i)) == 'specify') then
      if (moment(i) < min_moment(i)) then
        write(text,'(F8.4)') moment(i)
        text = adjustl(text)
        write(*,*) "Moment: "//trim(text)
        write(text,*) i
        text = adjustl(text)
        call ask_stop("Seed airfoil violates min_moment constraint for "//&
                      "operating point "//trim(text)//".")
      end if
    end if
  end do

! Evaluate objectives to establish scale factors for each point

  do i = 1, noppoint
    write(text,*) i
    text = adjustl(text)

    if (lift(i) <= 0.d0 .and. (trim(optimization_type(i)) == 'min-sink' .or.   &
        trim(optimization_type(i)) == 'max-glide') ) then
      write(*,*) "Error: operating point "//trim(text)//" has Cl <= 0. "//     &
                 "Cannot use "//trim(optimization_type(i))//" optimization "// &
                 "in this case."
      write(*,*)
      stop
    end if

    if (trim(optimization_type(i)) == 'min-sink') then
      checkval = drag(i)/lift(i)**1.5d0
    elseif (trim(optimization_type(i)) == 'max-glide') then
      checkval = drag(i)/lift(i)
    elseif (trim(optimization_type(i)) == 'min-drag') then
      checkval = drag(i)
    elseif (trim(optimization_type(i)) == 'max-lift') then
      checkval = 1.d0/lift(i)
    elseif (trim(optimization_type(i)) == 'max-xtr') then
      checkval = 1.d0/(0.5d0*(xtrt(i)+xtrb(i))+0.1d0)  ! Ensure no division by 0
    elseif (trim(optimization_type(i)) == 'max-lift-slope') then

!     Maximize dCl/dalpha (0.1 factor to ensure no division by 0)

      checkval = 0.d0
      if (i < noppoint) then
        if (alpha(i+1) > alpha(i)) then
          checkval = derv1f1(lift(i+1), lift(i),                               &
                             (alpha(i+1)-alpha(i)+0.1d0)*pi/180.d0)
        else
          checkval = derv1b1(lift(i+1), lift(i),                               &
                             (alpha(i)-alpha(i+1)+0.1d0)*pi/180.d0)
        end if
      end if

      if (i > 1) then
        if (alpha(i) > alpha(i-1)) then
          checkval = checkval + derv1b1(lift(i-1), lift(i),                    &
                                        (alpha(i)-alpha(i-1)+0.1d0)*pi/180.d0) 
        else
          checkval = checkval + derv1f1(lift(i-1), lift(i),                    &
                                        (alpha(i-1)-alpha(i)+0.1d0)*pi/180.d0) 
        end if
      end if
      if ( (i < noppoint) .and. (i > 1) ) checkval = checkval/2.d0 

!     4*pi factor is to move singularity location. Without it, negative
!     objective function occurs for Cl_a < 0. With this factor, it occurs for
!     Cl_a < -4*pi (hopefully never encountered).

      checkval = 1.d0/(checkval + 4.d0*pi)
      
    else
      write(*,*)
      write(*,*) "Error: requested optimization_type for operating point "//   &
                 trim(text)//" not recognized."
      stop
    end if
    scale_factor(i) = 1.d0/checkval
  end do

end subroutine check_seed

!=============================================================================80
!
! Asks user to stop or continue
!
!=============================================================================80
subroutine ask_stop(message)

  character(*), intent(in) :: message

  character :: choice
  logical :: valid_choice

! Get user input

  valid_choice = .false.
  do while (.not. valid_choice)
  
    write(*,'(A)') 'Warning: '//trim(message)
    write(*,'(A)', advance='no') 'Continue anyway? (y/n): '
    read(*,'(A)') choice

    if ( (choice == 'y') .or. (choice == 'Y') ) then
      valid_choice = .true.
      choice = 'y'
    else if ( ( choice == 'n') .or. (choice == 'N') ) then
      valid_choice = .true.
      choice = 'n'
    else
      write(*,'(A)') 'Please enter y or n.'
      valid_choice = .false.
    end if

  end do

! Stop or continue

  write(*,*)
  if (choice == 'n') stop

end subroutine ask_stop
  
end module input_sanity
