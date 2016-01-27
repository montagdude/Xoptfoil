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

module airfoil_evaluation

! Sets up and evaluates the objective function for an airfoil design

  use vardef
  use xfoil_driver, only : xfoil_options_type, xfoil_geom_options_type

  implicit none

  public
  private :: aero_objective_function, matchfoil_objective_function

  type(xfoil_options_type) :: xfoil_options
  type(xfoil_geom_options_type) :: xfoil_geom_options

! Variables used to check that XFoil results are repeatable when needed

  double precision :: checktol = 0.2d0
  double precision, dimension(max_op_points) :: maxlift = -100.d0
  double precision, dimension(max_op_points) :: mindrag = 100.d0

  contains

!=============================================================================80
!
! Generic objective function.  Selects either aero_objective_function or
! matchfoil_objective_function depending on whether match_foils = .true. or
! not.
!
!=============================================================================80
function objective_function(designvars)

  double precision, dimension(:), intent(in) :: designvars
  double precision :: objective_function

  if (match_foils) then
    objective_function = matchfoil_objective_function(designvars)
  else
    objective_function = aero_objective_function(designvars)
  end if

end function objective_function

!=============================================================================80
!
!  Objective function
!
!  Input: design variables (modes for top and bottom shape functions)
!  Output: objective function value based on airfoil performance
!
!=============================================================================80
function aero_objective_function(designvars)

  use math_deps,          only : interp_vector, curvature
  use parameterization,   only : top_shape_function, bot_shape_function,       &
                                 create_airfoil
  use xfoil_driver,       only : run_xfoil
  use xfoil_inc,          only : AMAX

  double precision, dimension(:), intent(in) :: designvars
  double precision :: aero_objective_function

  double precision, dimension(max(size(xseedt,1),size(xseedb,1))) :: x_interp, &
                                                 zt_interp, zb_interp, thickness
  double precision, dimension(size(xseedt,1)) :: zt_new
  double precision, dimension(size(xseedb,1)) :: zb_new
  double precision, dimension(size(xseedt,1)+size(xseedb,1)-1) :: curv
  integer :: nmodest, nmodesb, nptt, nptb, i, dvtbnd1, dvtbnd2, dvbbnd1,       &
             dvbbnd2, ncheckpt, nptint
  double precision :: penaltyval
  double precision :: tegap, growth1, growth2, maxgrowth, len1, len2
  double precision :: panang1, panang2, maxpanang, heightfactor
  integer, dimension(noppoint) :: checkpt_list
  character(7), dimension(noppoint) :: opm_check
  double precision, dimension(noppoint) :: opp_check, re_check, ma_check 
  double precision, dimension(noppoint) :: fd_check
  double precision, dimension(noppoint) :: lift, drag, moment, viscrms
  double precision, dimension(noppoint) :: clcheck, cdcheck, cmcheck, rmscheck
  double precision, dimension(noppoint) :: actual_flap_degrees
  logical, dimension(noppoint) :: checkpt
  double precision :: increment, curv1, curv2
  integer :: nreversals, ndvs
  double precision :: gapallow, maxthick, ffact
  integer :: check_idx, flap_idx, dvcounter
  double precision, parameter :: eps = 1.0D-08

  nmodest = size(top_shape_function,1)
  nmodesb = size(bot_shape_function,1)
  nptt = size(xseedt,1)
  nptb = size(xseedb,1)

! Set modes for top and bottom surfaces

  dvtbnd1 = 1
  if (trim(shape_functions) == 'naca') then
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
                      designvars(dvtbnd1:dvtbnd2), designvars(dvbbnd1:dvbbnd2),&
                      zt_new, zb_new, shape_functions, symmetrical)

! Format coordinates in a single loop in derived type

  do i = 1, nptt
    curr_foil%x(i) = xseedt(nptt-i+1)
    curr_foil%z(i) = zt_new(nptt-i+1)
  end do
  do i = 1, nptb-1
    curr_foil%x(i+nptt) = xseedb(i+1)
    curr_foil%z(i+nptt) = zb_new(i+1)
  end do

! Check geometry before running Xfoil: growth rates, LE and TE angles, etc.

  penaltyval = 0.d0
  maxgrowth = 0.d0

  len1 = sqrt((curr_foil%x(2)-curr_foil%x(1))**2.d0 +                          &
              (curr_foil%z(2)-curr_foil%z(1))**2.d0)
  do i = 2, nptt + nptb - 2
    len2 = sqrt((curr_foil%x(i+1)-curr_foil%x(i))**2.d0 +                      &
                (curr_foil%z(i+1)-curr_foil%z(i))**2.d0)
    growth1 = len2/len1
    growth2 = len1/len2
    if (max(growth1,growth2) > maxgrowth) maxgrowth = max(growth1,growth2)
    len1 = len2
  end do

! Penalty for too large growth rate

  penaltyval = penaltyval + max(0.d0,maxgrowth-growth_allowed)/1.d0

! Penalty for too blunt leading edge

  panang1 = atan((zt_new(2)-zt_new(1))/(xseedt(2)-xseedt(1))) *                &
            180.d0/acos(-1.d0)
  panang2 = atan((zb_new(1)-zb_new(2))/(xseedb(2)-xseedb(1))) *                &
            180.d0/acos(-1.d0)
  maxpanang = max(panang2,panang1)
  penaltyval = penaltyval + max(0.d0,maxpanang-89.99d0)/0.01d0

! Penalty for too sharp leading edge

  penaltyval = penaltyval + max(0.d0,abs(panang1-panang2)-20.d0)/5.d0

! Interpolate bottom surface to xseedt points (to check thickness)

  if (xseedt(nptt) <= xseedb(nptb)) then
    nptint = nptt
    call interp_vector(xseedb, zb_new, xseedt, zb_interp(1:nptt))
    x_interp(1:nptt) = xseedt
    zt_interp(1:nptt) = zt_new  
  else
    nptint = nptb
    call interp_vector(xseedt, zt_new, xseedb, zt_interp(1:nptb))
    x_interp(1:nptb) = xseedb
    zb_interp(1:nptb) = zb_new
  end if

! Compute thickness parameters

  tegap = zt_new(nptt) - zb_new(nptb)
  maxthick = 0.d0
  heightfactor = tan(min_te_angle*acos(-1.d0)/180.d0/2.d0)

  do i = 2, nptint - 1

!   Thickness array and max thickness

    thickness(i) = zt_interp(i) - zb_interp(i)
    if (thickness(i) > maxthick) maxthick = thickness(i)

!   Check if thinner than specified wedge angle on back half of airfoil

    if (xseedt(i) > 0.5d0) then
      gapallow = tegap + 2.d0 * heightfactor * (x_interp(nptint) -             &
                                                x_interp(i))
      if (thickness(i) < gapallow)                                             &
        penaltyval = penaltyval + (gapallow - thickness(i))/0.001d0
    end if

  end do

! Penalties for max thickness too low or high

  penaltyval = penaltyval + max(0.d0,min_thickness-maxthick)/0.1d0
  penaltyval = penaltyval + max(0.d0,maxthick-max_thickness)/0.1d0

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

    penaltyval = penaltyval + max(0.d0,dble(nreversals-max_curv_reverse))

  end if

! Check that number of flap optimize points are correct

  ndvs = size(designvars,1)
  if (nflap_optimize /= (ndvs - dvbbnd2)) then
    write(*,*) "Wrong number of design variables for flap deflections."
    write(*,*) "Please report this bug."
    stop
  end if

! Get actual flap angles based on design variables
! Also add a penalty for flap deflections outside the specified bounds

  ffact = initial_perturb/(max_flap_degrees - min_flap_degrees)
  actual_flap_degrees(1:noppoint) = flap_degrees(1:noppoint)
  dvcounter = dvbbnd2 + 1
  do i = 1, nflap_optimize
    flap_idx = flap_optimize_points(i)
    actual_flap_degrees(flap_idx) = designvars(dvcounter)/ffact
    penaltyval = penaltyval +                                                  &
                 max(0.d0,actual_flap_degrees(flap_idx)-max_flap_degrees)
    penaltyval = penaltyval +                                                  &
                 max(0.d0,min_flap_degrees-actual_flap_degrees(flap_idx))
    dvcounter = dvcounter + 1
  end do

! Exit if geometry and flap angles don't check out

  if (penaltyval > eps) then
    aero_objective_function = penaltyval*1.0D+06
    return
  end if

! Analyze airfoil at requested operating conditions with Xfoil

  call run_xfoil(curr_foil, xfoil_geom_options, op_point(1:noppoint),          &
                 op_mode(1:noppoint), reynolds(1:noppoint), mach(1:noppoint),  &
                 use_flap, x_flap, y_flap, actual_flap_degrees(1:noppoint),    &
                 xfoil_options, lift, drag, moment, viscrms)

! Determine if points need to be checked for xfoil consistency

  ncheckpt = 0
  checkpt_list(:) = 0
  checkpt(:) = .false.
  do i = 1, noppoint

!   Don't check the very first design

    if (maxlift(1) == -100.d0) exit

    if ((lift(i) > (1.d0 + checktol)*maxlift(i)) .or.                          &
        (drag(i) < (1.d0 - checktol)*mindrag(i))) then

      checkpt(i) = .true.
      ncheckpt = ncheckpt + 1
      checkpt_list(i) = ncheckpt
      opm_check(ncheckpt) = op_mode(i)
      opp_check(ncheckpt) = op_point(i)
      ma_check(ncheckpt) = mach(i)
      fd_check(ncheckpt) = actual_flap_degrees(i)

!     Perturb Reynolds number slightly to check that XFoil result is 
!     repeatable

      re_check(ncheckpt) = 0.997d0*reynolds(i)
 
    end if

  end do

! Analyze airfoil at perturbed operating points to check for repeatability

  anychecked: if (ncheckpt > 0) then

    call run_xfoil(curr_foil, xfoil_geom_options, opp_check(1:ncheckpt),       &
                   opm_check(1:ncheckpt), re_check(1:ncheckpt),                &
                   ma_check(1:ncheckpt), use_flap, x_flap, y_flap,             &
                   fd_check(1:ncheckpt), xfoil_options, clcheck, cdcheck,      &
                   cmcheck, rmscheck)

!   Keep the more conservative of the two runs

    do i = 1, noppoint

      ischecked: if (checkpt(i)) then

        check_idx = checkpt_list(i)

        checklift: if (clcheck(check_idx) < lift(i)) then
          lift(i) = clcheck(check_idx)
        end if checklift

        checkdrag: if (cdcheck(check_idx) > drag(i)) then
          drag(i) = cdcheck(check_idx)
        end if checkdrag

        checkmoment: if (cmcheck(check_idx) < moment(i)) then
          moment(i) = cmcheck(check_idx)
        end if checkmoment

        checkrms: if (rmscheck(check_idx) > viscrms(i)) then
          viscrms(i) = rmscheck(check_idx)
        end if checkrms

      end if ischecked

    end do

  end if anychecked

! Get objective function contribution from aerodynamics (aero performance
! times normalized weight)

  aero_objective_function = 0.d0

  do i = 1, noppoint

!   Extra checks for really bad designs

    if (viscrms(i) >= 1.d0) then
      lift(i) = -0.1d0
      drag(i) = 1000.d0
      moment(i) = -10.d0
    end if

!   Objective function evaluation

    if (trim(optimization_type(i)) == 'min-sink') then

!     Maximize Cl^1.5/Cd

      if (lift(i) > 0.d0) then
        increment = drag(i)/lift(i)**1.5d0*scale_factor(i)
      else
        increment = 1.D9   ! Big penalty for lift <= 0
      end if

    elseif (trim(optimization_type(i)) == 'max-glide') then

!     Maximize Cl/Cd

      if (lift(i) > 0.d0) then
        increment = drag(i)/lift(i)*scale_factor(i)
      else
        increment = 1.D9   ! Big penalty for lift <= 0
      end if

    elseif (trim(optimization_type(i)) == 'min-drag') then

!     Minimize Cd

      increment = drag(i)*scale_factor(i)

    elseif (trim(optimization_type(i)) == 'max-lift') then

!     Maximize Cl (at given angle of attack)

      if (lift(i) > 0.d0) then
        increment = scale_factor(i)/lift(i)
      else
        increment = 1.D9   ! Big penalty for lift <= 0
      end if

    else

      write(*,*)
      write(*,*) "Error: requested optimization_type not recognized."
      stop

    end if

!   Add contribution to the objective function

    aero_objective_function = aero_objective_function + weighting(i)*increment

  end do

! Add penalty for unconverged points

  do i = 1, noppoint
    penaltyval = penaltyval + max(0.d0,viscrms(i)-1.0D-04)/1.0D-04
  end do

! Add penalty for too low moment

  do i = 1, noppoint
    if (trim(moment_constraint_type(i)) /= 'none') then
      penaltyval = penaltyval + max(0.d0,min_moment(i)-moment(i))/0.1d0
    end if
  end do

! Add penalty for too large panel angles

  maxpanang = AMAX
  penaltyval = penaltyval + max(0.0d0,maxpanang-25.d0)/5.d0

! Add all penalties to objective function, and make them very large

  aero_objective_function = aero_objective_function + penaltyval*1.0D+06

! Update maxlift and mindrag only if it is a good design

  if (penaltyval <= eps) then
    do i = 1, noppoint
!$omp critical
      if (lift(i) > maxlift(i)) maxlift(i) = lift(i)
      if (drag(i) < mindrag(i)) mindrag(i) = drag(i)
!$omp end critical
    end do
  end if

!Bug check
!if (aero_objective_function < 0.5) then
!  print *, "penaltyval: ", penaltyval
!  do i = 1, noppoint
!    print *, drag(i), checkpt(i)
!  end do
!end if

end function aero_objective_function

!=============================================================================80
!
! Objective function for matching one airfoil to another (for testing shape
! functions, optimization algorithms, etc.).  Assumes x-values of points line
! up; this should be handled before optimizing.
!
!=============================================================================80
function matchfoil_objective_function(designvars)

  use parameterization,   only : top_shape_function, bot_shape_function,       &
                                 create_airfoil
  use math_deps,          only : norm_2

  double precision, dimension(:), intent(in) :: designvars
  double precision :: matchfoil_objective_function

  double precision, dimension(size(xseedt,1)) :: zt_new
  double precision, dimension(size(xseedb,1)) :: zb_new
  integer :: nmodest, nmodesb, nptt, nptb, dvtbnd, dvbbnd

  nmodest = size(top_shape_function,1)
  nmodesb = size(bot_shape_function,1)
  nptt = size(xseedt,1)
  nptb = size(xseedb,1)

! Set modes for top and bottom surfaces

  if (trim(shape_functions) == 'naca') then
    dvtbnd = nmodest
    dvbbnd = nmodest + nmodesb
  else
    dvtbnd = nmodest*3
    dvbbnd = nmodest*3 + nmodesb*3
  end if

! Create top and bottom surfaces by perturbation of seed airfoil

  call create_airfoil(xseedt, zseedt, xseedb, zseedb, designvars(1:dvtbnd),    &
                      designvars(dvtbnd+1:dvbbnd), zt_new, zb_new,             &
                      shape_functions, .false.)

! Evaluate the new airfoil, not counting fixed LE and TE points

  matchfoil_objective_function = norm_2(zt_new(2:nptt-1) - zmatcht(2:nptt-1))
  matchfoil_objective_function = matchfoil_objective_function +                &
                                 norm_2(zb_new(2:nptb-1) - zmatchb(2:nptb-1))

end function matchfoil_objective_function

!=============================================================================80
!
! Generic function to write design for the seed airfoil. Used due to desire
! to keep converterfunc in optimization.F90 generic (i.e., only takes
! designvars as input).
!
!=============================================================================80
function write_function_seed(designvars)

  double precision, dimension(:), intent(in) :: designvars
  integer :: write_function_seed

  write_function_seed = write_function(designvars, .true.)

end function write_function_seed

!=============================================================================80
!
! Generic function to write designs for airfoils other than the seed. Used due
! to desire to keep converterfunc in optimization.F90 generic (i.e., only takes
! designvars as input).
!
!=============================================================================80
function write_function_nonseed(designvars)

  double precision, dimension(:), intent(in) :: designvars
  integer :: write_function_nonseed

  write_function_nonseed = write_function(designvars)

end function write_function_nonseed

!=============================================================================80
!
! Generic function to write designs. Selects either 
! write_airfoil_optimization_progress or write_matchfoil_optimization_progress
! depending on whether match_foils = .true. or not.
!
!=============================================================================80
function write_function(designvars, isseed)

  double precision, dimension(:), intent(in) :: designvars
  logical, optional, intent(in) :: isseed
  integer :: write_function

  if (match_foils) then
    write_function = write_matchfoil_optimization_progress(designvars, isseed)
  else
    write_function = write_airfoil_optimization_progress(designvars, isseed)
  end if

end function write_function

!=============================================================================80
!
! Writes airfoil coordinates and polars to files during optimization
!
!=============================================================================80
function write_airfoil_optimization_progress(designvars, isseed)

  use math_deps,          only : interp_vector 
  use parameterization,   only : top_shape_function, bot_shape_function,       &
                                 create_airfoil
  use xfoil_driver,       only : run_xfoil

  double precision, dimension(:), intent(in) :: designvars
  logical, optional, intent(in) :: isseed
  integer :: write_airfoil_optimization_progress

  double precision, dimension(size(xseedt,1)) :: zt_new
  double precision, dimension(size(xseedb,1)) :: zb_new
  integer :: nmodest, nmodesb, nptt, nptb, i, dvtbnd1, dvtbnd2, dvbbnd1,       &
             dvbbnd2 
  double precision, dimension(noppoint) :: lift, drag, moment, viscrms
  double precision, dimension(noppoint) :: actual_flap_degrees
  double precision :: ffact
  integer :: ndvs, flap_idx, dvcounter

  nmodest = size(top_shape_function,1)
  nmodesb = size(bot_shape_function,1)
  nptt = size(xseedt,1)
  nptb = size(xseedb,1)

! Set modes for top and bottom surfaces

  dvtbnd1 = 1
  if (trim(shape_functions) == 'naca') then
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
                      designvars(dvtbnd1:dvtbnd2), designvars(dvbbnd1:dvbbnd2),&
                      zt_new, zb_new, shape_functions, symmetrical)

! Format coordinates in a single loop in derived type

  do i = 1, nptt
    curr_foil%x(i) = xseedt(nptt-i+1)
    curr_foil%z(i) = zt_new(nptt-i+1)
  end do
  do i = 1, nptb-1
    curr_foil%x(i+nptt) = xseedb(i+1)
    curr_foil%z(i+nptt) = zb_new(i+1)
  end do

! Check that number of flap optimize points are correct

  ndvs = size(designvars,1)
  if (nflap_optimize /= (ndvs - dvbbnd2)) then
    write(*,*) "Wrong number of design variables for flap deflections."
    write(*,*) "Please report this bug."
    stop
  end if

! Get actual flap angles based on design variables

  ffact = initial_perturb/(max_flap_degrees - min_flap_degrees)
  actual_flap_degrees(1:noppoint) = flap_degrees(1:noppoint)
  dvcounter = dvbbnd2 + 1
  do i = 1, nflap_optimize
    flap_idx = flap_optimize_points(i)
    actual_flap_degrees(flap_idx) = designvars(dvcounter)/ffact
    dvcounter = dvcounter + 1
  end do

! Analyze airfoil at requested operating conditions with Xfoil

  call run_xfoil(curr_foil, xfoil_geom_options, op_point(1:noppoint),          &
                 op_mode(1:noppoint), reynolds(1:noppoint), mach(1:noppoint),  &
                 use_flap, x_flap, y_flap, actual_flap_degrees(1:noppoint),    &
                 xfoil_options, lift, drag, moment, viscrms)

! Write coordinates and polars to file

! Set return value (needed for compiler)

  write_airfoil_optimization_progress = 0

end function write_airfoil_optimization_progress

!=============================================================================80
!
! Writes airfoil coordinates to foil during optimization to match one airfoil
! to another.
!
!=============================================================================80
function write_matchfoil_optimization_progress(designvars, isseed)

  double precision, dimension(:), intent(in) :: designvars
  logical, optional, intent(in) :: isseed
  integer :: write_matchfoil_optimization_progress

! Set return value (needed for compiler)

  write_matchfoil_optimization_progress = 0

end function write_matchfoil_optimization_progress

end module airfoil_evaluation
