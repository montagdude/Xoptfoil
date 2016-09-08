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
! Objective function with option to not add penalty value (used for seed
! airfoil)
!
!=============================================================================80
function objective_function_nopenalty(designvars)

  double precision, dimension(:), intent(in) :: designvars
  double precision :: objective_function_nopenalty

  if (match_foils) then
    objective_function_nopenalty = matchfoil_objective_function(designvars)
  else
    objective_function_nopenalty =                                             &
                    aero_objective_function(designvars, include_penalty=.false.)
  end if

end function objective_function_nopenalty

!=============================================================================80
!
!  Objective function
!
!  Input: design variables (modes for top and bottom shape functions)
!         include_penalty: optional input to enable/disable penalty function
!  Output: objective function value based on airfoil performance
!
!=============================================================================80
function aero_objective_function(designvars, include_penalty)

  use math_deps,       only : interp_vector, curvature
  use parametrization, only : top_shape_function, bot_shape_function,          &
                              create_airfoil
  use xfoil_driver,    only : run_xfoil
  use xfoil_inc,       only : AMAX, CAMBR

  double precision, dimension(:), intent(in) :: designvars
  logical, intent(in), optional :: include_penalty
  double precision :: aero_objective_function

  double precision, dimension(max(size(xseedt,1),size(xseedb,1))) :: x_interp, &
                                               zt_interp, zb_interp, thickness
  double precision, dimension(size(xseedt,1)) :: zt_new
  double precision, dimension(size(xseedb,1)) :: zb_new
  double precision, dimension(size(xseedt,1)) :: curvt
  double precision, dimension(size(xseedb,1)) :: curvb
  integer :: nmodest, nmodesb, nptt, nptb, i, dvtbnd1, dvtbnd2, dvbbnd1,       &
             dvbbnd2, ncheckpt, nptint
  double precision :: penaltyval
  double precision :: tegap, growth1, growth2, maxgrowth, len1, len2
  double precision :: panang1, panang2, maxpanang, heightfactor
  integer, dimension(noppoint) :: checkpt_list
  character(7), dimension(noppoint) :: opm_check
  double precision, dimension(noppoint) :: opp_check, re_check, ma_check 
  double precision, dimension(noppoint) :: fd_check
  double precision, dimension(noppoint) :: lift, drag, moment, viscrms, alpha, &
                                           xtrt, xtrb
  double precision, dimension(noppoint) :: clcheck, cdcheck, cmcheck, rmscheck,&
                                           alcheck, xtrtcheck, xtrbcheck
  double precision, dimension(noppoint) :: actual_flap_degrees
  logical, dimension(noppoint) :: checkpt
  logical :: check
  double precision :: increment, curv1, curv2
  integer :: nreversalst, nreversalsb, ndvs
  double precision :: gapallow, maxthick, ffact
  integer :: check_idx, flap_idx, dvcounter
  double precision, parameter :: epsexit = 1.0D-04
  double precision, parameter :: epsupdate = 1.0D-08
  logical :: penalize

  nmodest = size(top_shape_function,1)
  nmodesb = size(bot_shape_function,1)
  nptt = size(xseedt,1)
  nptb = size(xseedb,1)

! Enable / disable penalty function

  penalize = .true.
  if (present(include_penalty)) then
    if (.not. include_penalty) penalize = .false.
  end if

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
      penaltyval = penaltyval + max(0.d0,gapallow-thickness(i))/0.001d0
    end if

  end do

! Penalties for max thickness too low or high

  penaltyval = penaltyval + max(0.d0,min_thickness-maxthick)/0.1d0
  penaltyval = penaltyval + max(0.d0,maxthick-max_thickness)/0.1d0

! Check for curvature reversals

  if (check_curvature) then

!   Compute curvature on top and bottom

    curvt = curvature(nptt, xseedt, zt_new)
    curvb = curvature(nptb, xseedb, zb_new)

!   Check number of reversals that exceed the threshold

    nreversalst = 0
    curv1 = 0.d0
    do i = 2, nptt - 1
      if (abs(curvt(i)) >= curv_threshold) then
        curv2 = curvt(i)
        if (curv2*curv1 < 0.d0) nreversalst = nreversalst + 1
        curv1 = curv2
      end if
    end do

    nreversalsb = 0
    curv1 = 0.d0
    do i = 2, nptb - 1
      if (abs(curvb(i)) >= curv_threshold) then
        curv2 = curvb(i)
        if (curv2*curv1 < 0.d0) nreversalsb = nreversalsb + 1
        curv1 = curv2
      end if
    end do

    penaltyval = penaltyval + max(0.d0,dble(nreversalst-max_curv_reverse_top))
    penaltyval = penaltyval + max(0.d0,dble(nreversalsb-max_curv_reverse_bot))

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

  if ( (penaltyval > epsexit) .and. penalize ) then
    aero_objective_function = penaltyval*1.0D+06
    return
  end if

! Analyze airfoil at requested operating conditions with Xfoil

  call run_xfoil(curr_foil, xfoil_geom_options, op_point(1:noppoint),          &
                 op_mode(1:noppoint), reynolds(1:noppoint), mach(1:noppoint),  &
                 use_flap, x_flap, y_flap, actual_flap_degrees(1:noppoint),    &
                 xfoil_options, lift, drag, moment, viscrms, alpha, xtrt, xtrb)

! Add penalty for too large panel angles

  penaltyval = penaltyval + max(0.0d0,AMAX-25.d0)/5.d0

! Add penalty for camber outside of constraints

  penaltyval = penaltyval + max(0.d0,CAMBR-max_camber)/0.025d0
  penaltyval = penaltyval + max(0.d0,min_camber-CAMBR)/0.025d0

! Exit if panel angles and camber constraints don't check out

  if ( (penaltyval > epsexit) .and. penalize ) then
    aero_objective_function = penaltyval*1.0D+06
    return
  end if

! Determine if points need to be checked for xfoil consistency

  ncheckpt = 0
  checkpt_list(:) = 0
  checkpt(:) = .false.
  do i = 1, noppoint

!   Don't check the very first design

    if (maxlift(1) == -100.d0) exit

!   Check when lift or drag values are suspect

    check = .false.
    if (trim(optimization_type(i)) == 'min-drag') then
      if (drag(i) < (1.d0 - checktol)*mindrag(i)) check = .true.
    else
      if ((lift(i) > (1.d0 + checktol)*maxlift(i)) .or.                        &
          (drag(i) < (1.d0 - checktol)*mindrag(i))) check = .true.
    end if

    if (check) then
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
                   cmcheck, rmscheck, alcheck, xtrtcheck, xtrbcheck)

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

        checkxtrt: if (xtrtcheck(check_idx) < xtrt(i)) then
          xtrt(i) = xtrtcheck(check_idx)
        end if checkxtrt

        checkxtrb: if (xtrbcheck(check_idx) < xtrb(i)) then
          xtrb(i) = xtrbcheck(check_idx)
        end if checkxtrb

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

    elseif (trim(optimization_type(i)) == 'max-xtr') then

!     Maximize laminar flow on top and bottom (0.1 factor to ensure no
!     division by 0)

      increment = scale_factor(i)/(0.5d0*(xtrt(i)+xtrb(i))+0.1d0)

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

! Add all penalties to objective function, and make them very large

  if (penalize) aero_objective_function =                                      &
                aero_objective_function + penaltyval*1.0D+06

! Update maxlift and mindrag only if it is a good design

  if (penaltyval <= epsupdate) then
    do i = 1, noppoint
!$omp critical
      if (lift(i) > maxlift(i)) maxlift(i) = lift(i)
      if (drag(i) < mindrag(i)) mindrag(i) = drag(i)
!$omp end critical
    end do
  end if

end function aero_objective_function

!=============================================================================80
!
! Objective function for matching one airfoil to another (for testing shape
! functions, optimization algorithms, etc.).  Assumes x-values of points line
! up; this should be handled before optimizing.
!
!=============================================================================80
function matchfoil_objective_function(designvars)

  use parametrization, only : top_shape_function, bot_shape_function,          &
                              create_airfoil
  use math_deps,       only : norm_2

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
! Generic function to write designs. Selects either 
! write_airfoil_optimization_progress or write_matchfoil_optimization_progress
! depending on whether match_foils = .true. or not.
!
!=============================================================================80
function write_function(designvars, designcounter)

  double precision, dimension(:), intent(in) :: designvars
  integer, intent(in) :: designcounter
  integer :: write_function

  if (match_foils) then
    write_function = write_matchfoil_optimization_progress(designvars,         &
                                                           designcounter)
  else
    write_function = write_airfoil_optimization_progress(designvars,           &
                                                         designcounter)
  end if

end function write_function

!=============================================================================80
!
! Writes airfoil coordinates and polars to files during optimization
!
!=============================================================================80
function write_airfoil_optimization_progress(designvars, designcounter)

  use math_deps,       only : interp_vector 
  use parametrization, only : top_shape_function, bot_shape_function,          &
                              create_airfoil
  use xfoil_driver,    only : run_xfoil, xfoil_geometry_info

  double precision, dimension(:), intent(in) :: designvars
  integer, intent(in) :: designcounter
  integer :: write_airfoil_optimization_progress

  double precision, dimension(size(xseedt,1)) :: zt_new
  double precision, dimension(size(xseedb,1)) :: zb_new
  integer :: nmodest, nmodesb, nptt, nptb, i, dvtbnd1, dvtbnd2, dvbbnd1,       &
             dvbbnd2 
  double precision, dimension(noppoint) :: alpha, lift, drag, moment, viscrms, &
                                           xtrt, xtrb
  double precision, dimension(noppoint) :: actual_flap_degrees
  double precision :: ffact, maxt, xmaxt, maxc, xmaxc
  integer :: ndvs, flap_idx, dvcounter
 
  character(100) :: foilfile, polarfile, text
  character(8) :: maxtchar, xmaxtchar, maxcchar, xmaxcchar
  integer :: foilunit, polarunit

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
                 xfoil_options, lift, drag, moment, viscrms, alpha, xtrt, xtrb)

! Get geometry info

  call xfoil_geometry_info(maxt, xmaxt, maxc, xmaxc)
  write(maxtchar,'(F8.5)') maxt
  maxtchar = adjustl(maxtchar)
  write(xmaxtchar,'(F8.5)') xmaxt
  xmaxtchar = adjustl(xmaxtchar)
  write(maxcchar,'(F8.5)') maxc
  maxcchar = adjustl(maxcchar)
  write(xmaxcchar,'(F8.5)') xmaxc
  xmaxcchar = adjustl(xmaxcchar)

! Set output file names and identifiers

  foilfile = trim(output_prefix)//'_design_coordinates.dat'
  polarfile = trim(output_prefix)//'_design_polars.dat'

  foilunit = 13
  polarunit = 14

! Open files and write headers, if necessary

  if (designcounter == 0) then

!   Header for coordinate file

    write(*,*) "Writing coordinates for seed airfoil to file "//               &
               trim(foilfile)//" ..."
    open(unit=foilunit, file=foilfile, status='replace')
    write(foilunit,'(A)') 'title="Airfoil coordinates"'
    write(foilunit,'(A)') 'variables="x" "z"'
    write(foilunit,'(A)') 'zone t="Seed airfoil, maxt='//trim(maxtchar)//&
                          ', xmaxt='//trim(xmaxtchar)//', maxc='//&
                          trim(maxcchar)//', xmaxc='//trim(xmaxcchar)//'"'

!   Header for polar file

    write(*,*) "Writing polars for seed airfoil to file "//                    &
               trim(polarfile)//" ..."
    open(unit=polarunit, file=polarfile, status='replace')
    write(polarunit,'(A)') 'title="Airfoil polars"'
    write(polarunit,'(A)') 'variables="alpha" "cl" "cd" "cm" "xtrt" "xtrb"'
    write(polarunit,'(A)') 'zone t="Seed airfoil polar"'

  else

!   Format design counter as string

    write(text,*) designcounter
    text = adjustl(text)

!   Open coordinate file and write zone header

    write(*,*) "  Writing coordinates for design number "//trim(text)//        &
               " to file "//trim(foilfile)//" ..."
    open(unit=foilunit, file=foilfile, status='old', position='append',        &
         err=900)
    write(foilunit,'(A)') 'zone t="Airfoil, maxt='//trim(maxtchar)//&
                          ', xmaxt='//trim(xmaxtchar)//', maxc='//&
                          trim(maxcchar)//', xmaxc='//trim(xmaxcchar)//'", '//&
                          'SOLUTIONTIME='//trim(text)

    ! Open polar file and write zone header
    
    write(*,*) "  Writing polars for design number "//trim(text)//             &
               " to file "//trim(polarfile)//" ..."
    open(unit=polarunit, file=polarfile, status='old', position='append',      &
         err=901)
    write(polarunit,'(A)') 'zone t="Polars", SOLUTIONTIME='//trim(text)

  end if

! Write coordinates to file

  do i = 1, nptt + nptb - 1
    write(foilunit,'(2F12.6)') curr_foil%x(i), curr_foil%z(i)
  end do

! Write polars to file

  do i = 1, noppoint
    write(polarunit,'(6ES14.6)') alpha(i), lift(i), drag(i), moment(i),        &
                                 xtrt(i), xtrb(i)
  end do

! Close output files

  close(foilunit)
  close(polarunit)

! Set return value (needed for compiler)

  write_airfoil_optimization_progress = 0
  return

! Warning if there was an error opening design_coordinates file

900 write(*,*) "Warning: unable to open "//trim(foilfile)//". Skipping ..."
  write_airfoil_optimization_progress = 1
  return

! Warning if there was an error opening design_coordinates file

901 write(*,*) "Warning: unable to open "//trim(polarfile)//". Skipping ..."
  write_airfoil_optimization_progress = 2
  return

end function write_airfoil_optimization_progress

!=============================================================================80
!
! Writes airfoil coordinates to foil during optimization to match one airfoil
! to another.
!
!=============================================================================80
function write_matchfoil_optimization_progress(designvars, designcounter)

  use parametrization, only : top_shape_function, bot_shape_function,          &
                              create_airfoil

  double precision, dimension(:), intent(in) :: designvars
  integer, intent(in) :: designcounter
  integer :: write_matchfoil_optimization_progress

  double precision, dimension(size(xseedt,1)) :: zt_new
  double precision, dimension(size(xseedb,1)) :: zb_new
  integer :: i, nmodest, nmodesb, nptt, nptb, dvtbnd, dvbbnd

  character(100) :: foilfile, text
  integer :: foilunit

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

! Format coordinates in a single loop in derived type

  do i = 1, nptt
    curr_foil%x(i) = xseedt(nptt-i+1)
    curr_foil%z(i) = zt_new(nptt-i+1)
  end do
  do i = 1, nptb-1
    curr_foil%x(i+nptt) = xseedb(i+1)
    curr_foil%z(i+nptt) = zb_new(i+1)
  end do

! Set output file names and identifiers

  foilfile = trim(output_prefix)//'_design_coordinates.dat'
  foilunit = 13

! Open file and write header, if necessary

  if (designcounter == 0) then

!   Header for coordinate file

    write(*,*) "Writing coordinates for seed airfoil to file "//               &
               trim(foilfile)//" ..."
    open(unit=foilunit, file=foilfile, status='replace')
    write(foilunit,'(A)') 'title="Airfoil coordinates"'
    write(foilunit,'(A)') 'variables="x" "z"'
    write(foilunit,'(A)') 'zone t="Seed airfoil"'

  else

!   Format designcounter as string

    write(text,*) designcounter
    text = adjustl(text)

!   Open coordinate file and write zone header

    write(*,*) "  Writing coordinates for design number "//trim(text)//        &
               " to file "//trim(foilfile)//" ..."
    open(unit=foilunit, file=foilfile, status='old', position='append',        &
         err=910)
    write(foilunit,'(A)') 'zone t="Airfoil", SOLUTIONTIME='//trim(text)

  end if

! Write coordinates to file

  do i = 1, nptt + nptb - 1
    write(foilunit,'(2F12.6)') curr_foil%x(i), curr_foil%z(i)
  end do

! Close output file

  close(foilunit)

! Set return value (needed for compiler)

  write_matchfoil_optimization_progress = 0
  return

! Warning if there was an error opening design_coordinates file

910 write(*,*) "Warning: unable to open "//trim(foilfile)//". Skipping ..."
  write_matchfoil_optimization_progress = 1
  return

end function write_matchfoil_optimization_progress

!=============================================================================80
!
! Cleans up unused designs written prior to a restart
!
!=============================================================================80
function write_function_restart_cleanup(restart_status, global_search,         &
                                        local_search)

  character(*), intent(in) :: restart_status, global_search, local_search
  integer :: write_function_restart_cleanup

  integer :: restunit, ioerr, step, designcounter, foilunit, polarunit,        &
             histunit, ncoord
  integer :: i, j
  double precision, dimension(:,:), allocatable :: x, z, alpha, lift, drag,    &
                                                   moment, xtrt, xtrb
  double precision, dimension(:), allocatable :: fmin, relfmin, rad
  character(150), dimension(:), allocatable :: zoneinfo
  character(100) :: restfile, foilfile, polarfile, text
  character(11) :: stepchar
  character(20) :: fminchar, radchar
  character(25) :: relfminchar

! Print status

  write(*,*) 'Cleaning up unused designs written after restart save ...'

  restunit = 12
  foilunit = 13
  polarunit = 14
  histunit = 15

! Read last written design from restart file

  if (trim(restart_status) == 'global_optimization') then
    if (trim(global_search) == 'particle_swarm') then
      restfile = 'restart_pso_'//trim(output_prefix)
    else if (trim(global_search) == 'genetic_algorithm') then
      restfile = 'restart_ga_'//trim(output_prefix)
    end if
  else
    if (trim(local_search) == 'simplex') then
      restfile = 'restart_simplex_'//trim(output_prefix)
    end if
  end if

  open(unit=restunit, file=restfile, status='old', form='unformatted',         &
       iostat=ioerr)
  if (ioerr /= 0) then
    write_function_restart_cleanup = 1
    return
  end if
  read(restunit) step
  read(restunit) designcounter
  close(restunit)

! Allocate size of data arrays

  ncoord = size(xseedt,1) + size(xseedb,1) - 1
  allocate(x(ncoord,designcounter+1))
  allocate(z(ncoord,designcounter+1))
  allocate(alpha(noppoint,designcounter+1))
  allocate(lift(noppoint,designcounter+1))
  allocate(drag(noppoint,designcounter+1))
  allocate(moment(noppoint,designcounter+1))
  allocate(xtrt(noppoint,designcounter+1))
  allocate(xtrb(noppoint,designcounter+1))
  allocate(zoneinfo(designcounter+1))
  allocate(fmin(step))
  allocate(relfmin(step))
  allocate(rad(step))

! Open coordinates file

  foilfile = trim(output_prefix)//'_design_coordinates.dat'
  open(unit=foilunit, file=foilfile, status='old', iostat=ioerr)
  if (ioerr /= 0) then
    write_function_restart_cleanup = 2
    return
  end if

! Skip file header

  read(foilunit,*)
  read(foilunit,*)

! Read coordinates for each airfoil

  do i = 1, designcounter + 1
  
!   Read zone header

    read(foilunit,'(A)') zoneinfo(i)

!   Read coordinates

    do j = 1, ncoord
      read(foilunit,'(2F14.6)') x(j,i), z(j,i)
    end do

  end do

! Close coordinates file

  close(foilunit)

! Re-write coordinates file without the unused designs

  open(unit=foilunit, file=foilfile, status='replace')
  write(foilunit,'(A)') 'title="Airfoil coordinates"'
  write(foilunit,'(A)') 'variables="x" "z"'
  do i = 0, designcounter
!   Write zone header

    write(foilunit,'(A)') trim(zoneinfo(i+1))

!   Write coordinates

    do j = 1, ncoord
      write(foilunit,'(2F12.6)') x(j,i+1), z(j,i+1)
    end do
  end do

! Close coordinates file

  close(foilunit)

! Open history file

  open(unit=histunit, file='optimization_history.dat', status='old',           &
       iostat=ioerr)
  if (ioerr /= 0) then
    write_function_restart_cleanup = 3
    return
  end if

! Skip file header

  read(histunit,*)

! Read optimizer data at each iteration

  do i = 1, step
    read(histunit,*) j, fmin(i), relfmin(i), rad(i)
  end do

! Close history file

  close(histunit)

! Re-write history file without the unused iterations

  open(unit=histunit, file='optimization_history.dat', status='replace')
  write(histunit,'(A)') "Iteration  Objective function  "//&
                        "% Improvement over seed  Design radius"
  do i = 1, step
    write(stepchar,'(I11)') i
    write(fminchar,'(F14.10)') fmin(i)
    write(relfminchar,'(F14.10)') relfmin(i)
    write(radchar,'(ES14.6)') rad(i)
    write(histunit,'(A11,A20,A25,A20)') adjustl(stepchar), adjustl(fminchar),  &
                                        adjustl(relfminchar), adjustl(radchar)
  end do

! Close history file

  close(histunit)

! Return now if we're matching airfoils (no aero data)

  if (match_foils) then
    deallocate(x)
    deallocate(z)
    deallocate(alpha)
    deallocate(lift)
    deallocate(drag)
    deallocate(moment)
    deallocate(xtrt)
    deallocate(xtrb)
    deallocate(fmin)
    deallocate(relfmin)
    deallocate(rad)
    write(*,*) 'Finished cleaning up unused designs.'
    write_function_restart_cleanup = 0
    return
  end if

! Open polars file

  polarfile = trim(output_prefix)//'_design_polars.dat'
  open(unit=polarunit, file=polarfile, status='old', iostat=ioerr)
  if (ioerr /= 0) then
    write_function_restart_cleanup = 4
    return
  end if

! Skip file header

  read(polarunit,*)
  read(polarunit,*)

! Read polars for each airfoil

  do i = 1, designcounter + 1
  
!   Skip zone header

    read(polarunit,*)

!   Read polars

    do j = 1, noppoint
      read(polarunit,'(6ES14.6)') alpha(j,i), lift(j,i), drag(j,i),            &
                                  moment(j,i), xtrt(j,i), xtrb(j,i)
    end do

  end do

! Close polars file

  close(polarunit)

! Re-write polars file without the unused designs

  open(unit=polarunit, file=polarfile, status='replace')
  write(polarunit,'(A)') 'title="Airfoil polars"'
  write(polarunit,'(A)') 'variables="cl" "cd"'
  do i = 0, designcounter
!   Write zone header

    if (i == 0) then
      write(polarunit,'(A)') 'zone t="Seed airfoil polar"'
    else
      write(text,*) i
      text = adjustl(text)
      write(polarunit,'(A)') 'zone t="Polars", SOLUTIONTIME='//trim(text)
    end if

!   Write polars

    do j = 1, noppoint
      write(polarunit,'(6ES14.6)') alpha(j,i+1), lift(j,i+1), drag(j,i+1),     &
                                   moment(j,i+1), xtrt(j,i+1), xtrb(j,i+1)
    end do
  end do

! Close polars file

  close(polarunit)

! Deallocate data arrays

  deallocate(x)
  deallocate(z)
  deallocate(lift)
  deallocate(drag)
  deallocate(moment)
  deallocate(zoneinfo)
  deallocate(xtrt)
  deallocate(xtrb)
  deallocate(fmin)
  deallocate(relfmin)
  deallocate(rad)

! Print status

  write(*,*) 'Finished cleaning up unused designs.'
  write(*,*)

  write_function_restart_cleanup = 0

end function write_function_restart_cleanup

end module airfoil_evaluation
