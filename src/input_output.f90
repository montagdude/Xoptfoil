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

module input_output

! Module with subroutines for reading and writing of files

  implicit none

  contains

!=============================================================================80
!
! Subroutine to read inputs from namelist file
!
!=============================================================================80
subroutine read_inputs(input_file, search_type, global_search, local_search,   &
                       seed_airfoil, airfoil_file, naca_digits, nfunctions_top,&
                       nfunctions_bot, initial_perturb, pso_options,           &
                       ds_options, matchfoil_file)

  use vardef
  use optimization,       only : pso_options_type, ds_options_type
  use airfoil_operations, only : my_stop
  use airfoil_evaluation, only : xfoil_options, xfoil_geom_options
 
  character(*), intent(in) :: input_file
  character(80), intent(out) :: search_type, global_search, local_search,      &
                                seed_airfoil, airfoil_file, matchfoil_file
  character(4), intent(out) :: naca_digits
  integer, intent(out) :: nfunctions_top, nfunctions_bot
  double precision, intent(out) :: initial_perturb
  type(pso_options_type), intent(out) :: pso_options
  type(ds_options_type), intent(out) :: ds_options

  logical :: viscous_mode, silent_mode, fix_unconverged, pso_feasible_init,    &
             reinitialize, write_designs
  integer :: pop, pso_nstop, pso_maxit, simplex_maxit, bl_maxit,               &
             npan, pso_feasible_init_attempts
  double precision :: pso_tol, simplex_tol, ncrit, xtript, xtripb, vaccel
  double precision :: cvpar, cterat, ctrrat, xsref1, xsref2, xpref1, xpref2
  double precision :: pso_feasible_limit
  integer :: i, iunit, ioerr
  character(30) :: text

  namelist /optimization_options/ search_type, global_search, local_search,    &
            seed_airfoil, airfoil_file, naca_digits, shape_functions,          &
            nfunctions_top, nfunctions_bot, initial_perturb, write_designs
  namelist /operating_conditions/ noppoint, op_mode, op_point, reynolds, mach, &
            weighting, optimization_type 
  namelist /constraints/ seed_violation_handling, min_thickness, max_thickness,&
                         moment_constraint_type, min_moment, min_te_angle,     &
                         check_curvature, max_curv_reverse, curv_threshold,    &
                         symmetrical
  namelist /particle_swarm_options/ pop, pso_tol, pso_nstop, pso_maxit,        &
                                    pso_feasible_init, pso_feasible_limit,     &
                                    pso_feasible_init_attempts
  namelist /simplex_options/ simplex_tol, simplex_maxit
  namelist /xfoil_run_options/ ncrit, xtript, xtripb, viscous_mode,            &
            silent_mode, bl_maxit, vaccel, fix_unconverged, reinitialize
  namelist /xfoil_paneling_options/ npan, cvpar, cterat, ctrrat, xsref1,       &
            xsref2, xpref1, xpref2
  namelist /matchfoil_options/ match_foils, matchfoil_file

! Open input file

  iunit = 12
  open(unit=iunit, file=input_file, status='old', iostat=ioerr)
  if (ioerr /= 0) then
    write(*,*)
    write(*,*) 'Error: could not find input file '//trim(input_file)//'.'
    write(*,*)
    stop
  end if

! Set defaults for main namelist options

  search_type = 'global_and_local'
  global_search = 'particle_swarm'
  local_search = 'simplex'
  seed_airfoil = 'four_digit'
  naca_digits = '0012'
  shape_functions = 'naca'
  nfunctions_top = 10
  nfunctions_bot = 10
  initial_perturb = 0.025d0
  write_designs = .true.

! Read main namelist options

  read(iunit, nml=optimization_options)

! Set defaults for operating conditions and constraints

  noppoint = 1
  op_mode(:) = 'spec-cl'
  op_point(:) = 0.d0
  reynolds(:) = 1.0D+05
  mach(:) = 0.0
  weighting(:) = 1.d0
  optimization_type(:) = 'min-drag'

  seed_violation_handling = 'stop'
  min_thickness = 0.06d0
  max_thickness = 1000.d0
  moment_constraint_type = 'use_seed'
  min_moment = -1.d0
  min_te_angle = 5.d0
  check_curvature = .false.
  max_curv_reverse = 3
  curv_threshold = 0.30d0
  symmetrical = .false.

! Read operating conditions and constraints

  read(iunit, nml=operating_conditions)
  read(iunit, nml=constraints)

! Normalize weightings for operating points

  weighting = weighting/sum(weighting(1:noppoint))

! Error checking and setting search algorithm options

  if (trim(search_type) /= 'global_and_local' .and. trim(search_type) /=       &
      'global' .and. trim(search_type) /= 'local') then
    write(*,*)
    write(*,*) "Error: search_type must be 'global_and_local', 'global', "//   &
               "or 'local'."
    write(*,*)
    stop
  end if

! Set default particle swarm options

  pop = 40
  pso_tol = 1.0D-04
  pso_nstop = 40
  pso_maxit = 300
  pso_feasible_init = .true.
  pso_feasible_limit = 5.0D+04
  pso_feasible_init_attempts = 100

! Set default simplex search options

  simplex_tol = 1.0D-05
  simplex_maxit = 1000

  if (trim(search_type) == 'global_and_local' .or. trim(search_type) ==        &
      'global') then
  
    if(trim(global_search) == 'particle_swarm') then

!     Set default PSO options

      if (trim(shape_functions) == 'naca') then
        pso_options%const = .false.
      else
        pso_options%const = .true.
      end if

!     Read PSO options and put them into derived type

      read(iunit, nml=particle_swarm_options)
      pso_options%pop = pop
      pso_options%tol = pso_tol
      pso_options%nstop = pso_nstop
      pso_options%maxspeed = initial_perturb
      pso_options%maxit = pso_maxit
      pso_options%feasible_init = pso_feasible_init
      pso_options%feasible_limit = pso_feasible_limit
      pso_options%feasible_init_attempts = pso_feasible_init_attempts
      pso_options%write_designs = write_designs

    else

      write(*,*)
      write(*,*) "Namelist error: global search type '"//trim(global_search)// &
                 "' is not available."
      write(*,*)
      stop
     
    end if

  end if

  if (trim(search_type) == 'global_and_local' .or. trim(search_type) ==        &
      'local') then

    if (trim(local_search) == 'simplex') then

!     Read simplex search options and put them into derived type

      read(iunit, nml=simplex_options)
      ds_options%tol = simplex_tol
      ds_options%maxit = simplex_maxit
      ds_options%write_designs = write_designs

    else

      write(*,*)
      write(*,*) "Namelist error: local search type '"//trim(local_search)//   &
                 "' is not available."
      write(*,*)
      stop
     
    end if

  end if 

! Set default xfoil aerodynamics and paneling options

  ncrit = 9.d0
  xtript = 1.d0
  xtripb = 1.d0
  viscous_mode = .true.
  silent_mode = .true.
  bl_maxit = 100
  vaccel = 0.01d0
  fix_unconverged = .true.
  reinitialize = .true.

  npan = 160
  cvpar = 1.d0
  cterat = 0.15d0
  ctrrat = 0.2d0
  xsref1 = 1.d0
  xsref2 = 1.d0
  xpref1 = 1.d0
  xpref2 = 1.d0

! Read xfoil options and put them into derived types

  read(iunit, nml=xfoil_run_options)
  read(iunit, nml=xfoil_paneling_options)

  xfoil_options%ncrit = ncrit
  xfoil_options%xtript = xtript
  xfoil_options%xtripb = xtripb
  xfoil_options%viscous_mode = viscous_mode
  xfoil_options%silent_mode = silent_mode
  xfoil_options%maxit = bl_maxit
  xfoil_options%vaccel = vaccel
  xfoil_options%fix_unconverged = fix_unconverged
  xfoil_options%reinitialize = reinitialize

  xfoil_geom_options%npan = npan
  xfoil_geom_options%cvpar = cvpar
  xfoil_geom_options%cterat = cterat
  xfoil_geom_options%ctrrat = ctrrat
  xfoil_geom_options%xsref1 = xsref1
  xfoil_geom_options%xsref2 = xsref2
  xfoil_geom_options%xpref1 = xpref1
  xfoil_geom_options%xpref2 = xpref2

! Option to match seed airfoil to another instead of aerodynamic optimization

  match_foils = .false.
  read(iunit, nml=matchfoil_options)

! Close the input file

  close(iunit)

! Echo namelist options for checking purposes

  write(*,*)
  write(*,*) 'Echoing program options:'
  write(*,*)

! Optimization options namelist

  write(*,'(A)') " &optimization_options"
  write(*,*) " search_type = '"//trim(search_type)//"'"
  write(*,*) " global_search = '"//trim(global_search)//"'"
  write(*,*) " local_search = '"//trim(local_search)//"'"
  write(*,*) " seed_airfoil = '"//trim(seed_airfoil)//"'"
  write(*,*) " airfoil_file = '"//trim(airfoil_file)//"'"
  write(*,*) " naca_digits = '"//trim(naca_digits)//"'"
  write(*,*) " shape_functions = '"//trim(shape_functions)//"'"
  write(*,*) " nfunctions_top = ", nfunctions_top
  write(*,*) " nfunctions_bot = ", nfunctions_bot
  write(*,*) " initial_perturb = ", initial_perturb
  write(*,*) " write_designs = ", write_designs
  write(*,'(A)') " /"
  write(*,*)

! Operating conditions namelist

  write(*,'(A)') " &operating_conditions"
  write(*,*) " noppoint = ", noppoint
  do i = 1, noppoint
    write(text,*) i
    text = adjustl(text)
    write(*,*) " optimization_type("//trim(text)//") = '"//                    &
               trim(optimization_type(i))//"'"
    write(*,*) " op_mode("//trim(text)//") = '"//trim(op_mode(i))//"'"
    write(*,*) " op_point("//trim(text)//") = ", op_point(i)
    write(*,'(A,es17.8)') "  reynolds("//trim(text)//") = ", reynolds(i)
    write(*,*) " mach("//trim(text)//") = ", mach(i)
    write(*,*) " weighting("//trim(text)//") = ", weighting(i)
  end do
  write(*,'(A)') " /"
  write(*,*)

! Constraints namelist

  write(*,'(A)') " &constraints"
  write(*,*) " seed_violation_handling = "//trim(seed_violation_handling)
  write(*,*) " min_thickness = ", min_thickness
  write(*,*) " max_thickness = ", max_thickness
  write(*,*) " moment_constraint_type = "//trim(moment_constraint_type)
  write(*,*) " min_moment = ", min_moment
  write(*,*) " min_te_angle = ", min_te_angle
  write(*,*) " check_curvature = ", check_curvature
  write(*,*) " max_curv_reverse = ", max_curv_reverse
  write(*,*) " curv_threshold = ", curv_threshold
  write(*,*) " symmetrical = ", symmetrical
  write(*,'(A)') " /"
  write(*,*)

  if (trim(search_type) == 'global_and_local' .or. trim(search_type) ==        &
      'global') then

    if(trim(global_search) == 'particle_swarm') then

!     Particle swarm namelist

      write(*,'(A)') " &particle_swarm_options"
      write(*,*) " pop = ", pso_options%pop
      write(*,*) " pso_tol = ", pso_options%tol
      write(*,*) " pso_nstop = ", pso_options%nstop
      write(*,*) " pso_maxit = ", pso_options%maxit
      write(*,*) " pso_feasible_init = ", pso_options%feasible_init
      write(*,*) " pso_feasible_limit = ", pso_options%feasible_limit
      write(*,*) " pso_feasible_init_attempts = ",                             &
                   pso_options%feasible_init_attempts
      write(*,'(A)') " /"
      write(*,*)

    end if

  end if

  if (trim(search_type) == 'global_and_local' .or. trim(search_type) ==        &
      'local') then

    if(trim(local_search) == 'simplex') then

!     Simplex search namelist

      write(*,'(A)') " &simplex_options"
      write(*,*) " simplex_tol = ", ds_options%tol
      write(*,*) " simplex_maxit = ", ds_options%maxit
      write(*,'(A)') " /"
      write(*,*)

    end if

  end if

! Xfoil run options namelist

  write(*,'(A)') " &xfoil_run_options"
  write(*,*) " ncrit = ", xfoil_options%ncrit
  write(*,*) " xtript = ", xfoil_options%xtript
  write(*,*) " xtripb = ", xfoil_options%xtripb
  write(*,*) " viscous_mode = ", xfoil_options%viscous_mode
  write(*,*) " silent_mode = ", xfoil_options%silent_mode
  write(*,*) " bl_maxit = ", xfoil_options%maxit
  write(*,*) " vaccel = ", xfoil_options%vaccel
  write(*,*) " fix_unconverged = ", xfoil_options%fix_unconverged
  write(*,*) " reinitialize = ", xfoil_options%reinitialize
  write(*,'(A)') " /"
  write(*,*)

! Xfoil paneling options namelist

  write(*,'(A)') " &xfoil_paneling_options"
  write(*,*) " npan = ", xfoil_geom_options%npan
  write(*,*) " cvpar = ", xfoil_geom_options%cvpar
  write(*,*) " cterat = ", xfoil_geom_options%cterat
  write(*,*) " ctrrat = ", xfoil_geom_options%ctrrat
  write(*,*) " xsref1 = ", xfoil_geom_options%xsref1
  write(*,*) " xsref2 = ", xfoil_geom_options%xsref2
  write(*,*) " xpref1 = ", xfoil_geom_options%xpref1
  write(*,*) " xpref2 = ", xfoil_geom_options%xpref2
  write(*,'(A)') " /"
  write(*,*)

! Matchfoil options

  write(*,'(A)') " &matchfoil_options"
  write(*,*) " match_foils = ", match_foils
  write(*,*) " matchfoil_file = '"//trim(matchfoil_file)//"'"
  write(*,'(A)') " /"
  write(*,*)

! Check that inputs are reasonable

! Optimization settings

  if (trim(seed_airfoil) /= 'from_file' .and.                                  &
      trim(seed_airfoil) /= 'four_digit')                                      &
    call my_stop("seed_airfoil must be 'from_file' or 'four_digit'.")

  if (trim(shape_functions) /= 'hicks-henne' .and.                             &
      trim(shape_functions) /= 'naca')                                         &
    call my_stop("shape_functions must be 'hicks-henne' or 'naca'.")

  if (nfunctions_top < 1)                                                      &
    call my_stop("nfunctions_top must be > 0.")

  if (nfunctions_bot < 1)                                                      &
    call my_stop("nfunctions_bot must be > 0.")

  if (initial_perturb <= 0.d0)                                                 &
    call my_stop("initial_perturb must be > 0.")

! Operating points

  if (noppoint < 1) call my_stop("noppoint must be > 0.")

  do i = 1, noppoint
    if (trim(op_mode(i)) /= 'spec-cl' .and. trim(op_mode(i)) /= 'spec-al')     &
      call my_stop("op_mode must be 'spec-al' or 'spec-cl'.")
    if (reynolds(i) <= 0.d0) call my_stop("reynolds must be > 0.")
    if (mach(i) < 0.d0) call my_stop("mach must be >= 0.")
    if (weighting(i) <= 0.d0) call my_stop("weighting must be > 0.")
    if (trim(optimization_type(i)) /= 'min-drag' .and.                         &
      trim(optimization_type(i)) /= 'max-glide' .and.                          &
      trim(optimization_type(i)) /= 'min-sink' .and.                           &
      trim(optimization_type(i)) /= 'max-lift')                                &
      call my_stop("optimization_type must be 'min-drag', 'max-glide', "//     &
                   "min-sink', or 'max-lift'.")
  end do

! Constraints

  if (trim(seed_violation_handling) /= 'stop' .and.                            &
      trim(seed_violation_handling) /= 'warn')                                 &
    call my_stop("seed_violation_handling must be 'stop' or 'warn'.")
  if (min_thickness <= 0.d0) call my_stop("min_thickness must be > 0.")
  if (max_thickness <= 0.d0) call my_stop("max_thickness must be > 0.")
  if (trim(moment_constraint_type) /= 'use_seed' .and.                         &
      trim(moment_constraint_type) /= 'specify' .and.                          &
      trim(moment_constraint_type) /= 'none')                                  &
    call my_stop("moment_constraint_type must be 'use_seed', 'specify', or "// &
                 "'none'.")
  if (min_te_angle <= 0.d0) call my_stop("min_te_angle must be > 0.")
  if (max_curv_reverse < 1) call my_stop("max_curv_reverse must be > 0.")
  if (curv_threshold <= 0.d0) call my_stop("curv_threshold must be > 0.")
  if (symmetrical)                                                             &
    write(*,*) "Mirroring top half of seed airfoil for symmetrical constraint."
    
! Particle swarm options

  if (pop < 1) call my_stop("pop must be > 0.")
  if (pso_tol <= 0.d0) call my_stop("pso_tol must be > 0.")
  if (pso_nstop < 1) call my_stop("pso_nstop must be > 0.")
  if (pso_maxit < 1) call my_stop("pso_maxit must be > 0.")  
  if ((pso_feasible_limit <= 0.d0) .and. pso_feasible_init)                    &
    call my_stop("pso_feasible_limit must be > 0.")
  if ((pso_feasible_init_attempts < 1) .and. pso_feasible_init)                &
    call my_stop("pso_feasible_init_attempts must be > 0.")

! Simplex options

  if (simplex_tol <= 0.d0) call my_stop("simplex_tol must be > 0.")
  if (simplex_maxit < 1) call my_stop("simplex_maxit must be > 0.")  

! XFoil run options

  if (ncrit < 0.d0) call my_stop("ncrit must be >= 0.")
  if (xtript < 0.d0 .or. xtript > 1.d0)                                        &
    call my_stop("xtript must be >= 0. and <= 1.")
  if (xtripb < 0.d0 .or. xtripb > 1.d0)                                        &
    call my_stop("xtripb must be >= 0. and <= 1.")
  if (bl_maxit < 1) call my_stop("bl_maxit must be > 0.")
  if (vaccel < 0.d0) call my_stop("vaccel must be >= 0.")
  
! XFoil paneling options

  if (npan < 20) call my_stop("npan must be >= 20.")
  if (cvpar <= 0.d0) call my_stop("cvpar must be > 0.")
  if (cterat <= 0.d0) call my_stop("cterat must be > 0.")
  if (ctrrat <= 0.d0) call my_stop("ctrrat must be > 0.")
  if (xsref1 < 0.d0) call my_stop("xsref1 must be >= 0.")
  if (xsref2 < xsref1) call my_stop("xsref2 must be >= xsref1")
  if (xsref2 > 1.d0) call my_stop("xsref2 must be <= 1.")
  if (xpref1 < 0.d0) call my_stop("xpref1 must be >= 0.")
  if (xpref2 < xpref1) call my_stop("xpref2 must be >= xpref1")
  if (xpref2 > 1.d0) call my_stop("xpref2 must be <= 1.")

end subroutine read_inputs

!=============================================================================80
!
! Subroutine to read inputs from namelist file - for xfoil_only
!
!=============================================================================80
subroutine read_inputs_xfoil_only(airfoil_file)

  use vardef,             only : max_op_points, noppoint, op_mode, op_point,   &
                                 reynolds, mach 
  use airfoil_evaluation, only : xfoil_options, xfoil_geom_options
 
  character(80), intent(out) :: airfoil_file

  logical :: viscous_mode, silent_mode, fix_unconverged, reinitialize
  integer :: bl_maxit, npan
  double precision :: ncrit, xtript, xtripb, vaccel
  double precision :: cvpar, cterat, ctrrat, xsref1, xsref2, xpref1, xpref2
  integer :: i, iunit, ioerr
  character(30) :: text

  namelist /airfoil_to_load/ airfoil_file
  namelist /operating_conditions/ noppoint, op_mode, op_point, reynolds, mach
  namelist /xfoil_run_options/ ncrit, xtript, xtripb, viscous_mode,            &
            silent_mode, bl_maxit, vaccel, fix_unconverged, reinitialize
  namelist /xfoil_paneling_options/ npan, cvpar, cterat, ctrrat, xsref1,       &
            xsref2, xpref1, xpref2

! Open input file

  iunit = 12
  open(unit=iunit, file='inputs_xfoil_only.txt', status='old', iostat=ioerr)
  if (ioerr /= 0) then
    write(*,*)
    write(*,*) 'Error: could not find input file inputs_xfoil_only.txt.'
    write(*,*)
    stop
  end if

! Read airfoil_to_load namelist options

  read(iunit, nml=airfoil_to_load)

! Set defaults for operating conditions

  noppoint = 1
  op_mode(:) = 'spec-cl'
  op_point(:) = 0.d0
  reynolds(:) = 1.0D+05
  mach(:) = 0.0

! Read operating conditions and constraints

  read(iunit, nml=operating_conditions)

! Set default xfoil aerodynamics and paneling options

  xfoil_options%ncrit = 9.d0
  xfoil_options%xtript = 1.d0
  xfoil_options%xtripb = 1.d0
  xfoil_options%viscous_mode = .true.
  xfoil_options%silent_mode = .true.
  xfoil_options%maxit = 100
  xfoil_options%vaccel = 0.01d0
  xfoil_options%fix_unconverged = .true.
  xfoil_options%reinitialize = .true.

  xfoil_geom_options%npan = 160
  xfoil_geom_options%cvpar = 1.d0
  xfoil_geom_options%cterat = 0.15d0
  xfoil_geom_options%ctrrat = 0.2d0
  xfoil_geom_options%xsref1 = 1.d0
  xfoil_geom_options%xsref2 = 1.d0
  xfoil_geom_options%xpref1 = 1.d0
  xfoil_geom_options%xpref2 = 1.d0

! Read xfoil options and put them into derived types

  read(iunit, nml=xfoil_run_options)
  read(iunit, nml=xfoil_paneling_options)

  xfoil_options%ncrit = ncrit
  xfoil_options%xtript = xtript
  xfoil_options%xtripb = xtripb
  xfoil_options%viscous_mode = viscous_mode
  xfoil_options%silent_mode = silent_mode
  xfoil_options%maxit = bl_maxit
  xfoil_options%vaccel = vaccel
  xfoil_options%fix_unconverged = fix_unconverged
  xfoil_options%reinitialize = reinitialize

  xfoil_geom_options%npan = npan
  xfoil_geom_options%cvpar = cvpar
  xfoil_geom_options%cterat = cterat
  xfoil_geom_options%ctrrat = ctrrat
  xfoil_geom_options%xsref1 = xsref1
  xfoil_geom_options%xsref2 = xsref2
  xfoil_geom_options%xpref1 = xpref1
  xfoil_geom_options%xpref2 = xpref2

! Close the input file

  close(iunit)

! Echo namelist options for checking purposes

  write(*,*)
  write(*,*) 'Echoing program options:'
  write(*,*)

! airfoil_to_load namelist

  write(*,'(A)') " &airfoil_to_load"
  write(*,*) " airfoil_file = '"//trim(airfoil_file)//"'"
  write(*,'(A)') " /"
  write(*,*)

! Operating conditions namelist

  write(*,'(A)') " &operating_conditions"
  write(*,*) " noppoint = ", noppoint
  do i = 1, noppoint
    write(text,*) i
    text = adjustl(text)
    write(*,*) " op_mode("//trim(text)//") = '"//trim(op_mode(i))//"'"
    write(*,*) " op_point("//trim(text)//") = ", op_point(i)
    write(*,'(A,es17.8)') "  reynolds("//trim(text)//") = ", reynolds(i)
    write(*,*) " mach("//trim(text)//") = ", mach(i)
  end do
  write(*,'(A)') " /"
  write(*,*)

! Xfoil run options namelist

  write(*,'(A)') " &xfoil_run_options"
  write(*,*) " ncrit = ", xfoil_options%ncrit
  write(*,*) " xtript = ", xfoil_options%xtript
  write(*,*) " xtripb = ", xfoil_options%xtripb
  write(*,*) " viscous_mode = ", xfoil_options%viscous_mode
  write(*,*) " silent_mode = ", xfoil_options%silent_mode
  write(*,*) " bl_maxit = ", xfoil_options%maxit
  write(*,*) " vaccel = ", xfoil_options%vaccel
  write(*,*) " fix_unconverged = ", xfoil_options%fix_unconverged
  write(*,*) " reinitialize = ", xfoil_options%reinitialize
  write(*,'(A)') " /"
  write(*,*)

! Xfoil paneling options namelist

  write(*,'(A)') " &xfoil_paneling_options"
  write(*,*) " npan = ", xfoil_geom_options%npan
  write(*,*) " cvpar = ", xfoil_geom_options%cvpar
  write(*,*) " cterat = ", xfoil_geom_options%cterat
  write(*,*) " ctrrat = ", xfoil_geom_options%ctrrat
  write(*,*) " xsref1 = ", xfoil_geom_options%xsref1
  write(*,*) " xsref2 = ", xfoil_geom_options%xsref2
  write(*,*) " xpref1 = ", xfoil_geom_options%xpref1
  write(*,*) " xpref2 = ", xfoil_geom_options%xpref2
  write(*,'(A)') " /"
  write(*,*)

end subroutine read_inputs_xfoil_only

!=============================================================================80
!
! Reads command line arguments for input file name and output file prefix
!
!=============================================================================80
subroutine read_clo(input_file, output_prefix)

  character(*), intent(inout) :: input_file, output_prefix

! Set default names

  input_file = "inputs.txt"
  output_prefix = "optfoil"

! Read supplied names

  if (iargc() >= 1) then
    call getarg(1, input_file)
  end if
  if (iargc() >= 2) then
    call getarg(2, output_prefix)
  end if

end subroutine read_clo

end module input_output
