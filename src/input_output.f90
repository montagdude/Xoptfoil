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
                       nfunctions_bot, restart, restart_write_freq,            &
                       constrained_dvs, pso_options, ga_options, ds_options,   &
                       matchfoil_file)

  use vardef
  use particle_swarm,     only : pso_options_type
  use genetic_algorithm,  only : ga_options_type
  use simplex_search,     only : ds_options_type
  use airfoil_operations, only : my_stop
  use airfoil_evaluation, only : xfoil_options, xfoil_geom_options
 
  character(*), intent(in) :: input_file
  character(80), intent(out) :: search_type, global_search, local_search,      &
                                seed_airfoil, airfoil_file, matchfoil_file
  character(4), intent(out) :: naca_digits
  integer, intent(out) :: nfunctions_top, nfunctions_bot
  integer, dimension(:), allocatable, intent(inout) :: constrained_dvs
  type(pso_options_type), intent(out) :: pso_options
  type(ga_options_type), intent(out) :: ga_options
  type(ds_options_type), intent(out) :: ds_options

  logical :: viscous_mode, silent_mode, fix_unconverged, feasible_init,        &
             reinitialize, restart, write_designs
  integer :: restart_write_freq, pso_pop, pso_maxit, simplex_maxit, bl_maxit,  &
             npan, feasible_init_attempts
  integer :: ga_pop, ga_maxit
  double precision :: pso_tol, simplex_tol, ncrit, xtript, xtripb, vaccel
  double precision :: cvpar, cterat, ctrrat, xsref1, xsref2, xpref1, xpref2
  double precision :: feasible_limit
  double precision :: ga_tol, parent_fraction, roulette_selection_pressure,    &
                      tournament_fraction, crossover_range_factor,             &
                      mutant_probability, chromosome_mutation_rate,            &
                      mutation_range_factor
  integer :: i, iunit, ioerr, iostat1, counter, idx
  character(30) :: text
  character(10) :: pso_convergence_profile, parents_selection_method
  character(200) :: msg

  namelist /optimization_options/ search_type, global_search, local_search,    &
            seed_airfoil, airfoil_file, naca_digits, shape_functions,          &
            nfunctions_top, nfunctions_bot, initial_perturb, min_bump_width,   &
            restart, restart_write_freq, write_designs
  namelist /operating_conditions/ noppoint, op_mode, op_point, reynolds, mach, &
            use_flap, x_flap, y_flap, flap_selection, flap_degrees, weighting, &
            optimization_type 
  namelist /constraints/ seed_violation_handling, min_thickness, max_thickness,&
                         moment_constraint_type, min_moment, min_te_angle,     &
                         check_curvature, max_curv_reverse, curv_threshold,    &
                         symmetrical, min_flap_degrees, max_flap_degrees
  namelist /initialization/ feasible_init, feasible_limit,                     &
                            feasible_init_attempts
  namelist /particle_swarm_options/ pso_pop, pso_tol, pso_maxit,               &
                                    pso_convergence_profile
  namelist /genetic_algorithm_options/ ga_pop, ga_tol, ga_maxit,               &
            parents_selection_method, parent_fraction,                         &
            roulette_selection_pressure, tournament_fraction,                  &
            crossover_range_factor, mutant_probability,                        &
            chromosome_mutation_rate, mutation_range_factor
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
  shape_functions = 'hicks-henne'
  min_bump_width = 0.5d0
  nfunctions_top = 4
  nfunctions_bot = 4
  initial_perturb = 0.025d0
  restart = .false.
  restart_write_freq = 20
  write_designs = .true.

! Read main namelist options

  rewind(iunit)
  read(iunit, iostat=iostat1, nml=optimization_options)
  if (iostat1 /= 0) then
    msg = "Unrecognized variable name in namelist "//&
    "'optimization_options.' Check User Guide for available variables."
    call my_stop(msg)
  end if

! Error checking and setting search algorithm options

  if (trim(search_type) /= 'global_and_local' .and. trim(search_type) /=       &
      'global' .and. trim(search_type) /= 'local') then
    write(*,*)
    write(*,*) "Error: search_type must be 'global_and_local', 'global', "//   &
               "or 'local'."
    write(*,*)
    stop
  end if

! Set defaults for operating conditions and constraints

  noppoint = 1
  use_flap = .false.
  x_flap = 0.75d0
  y_flap = 0.d0
  op_mode(:) = 'spec-cl'
  op_point(:) = 0.d0
  optimization_type(:) = 'min-drag'
  reynolds(:) = 1.0D+05
  mach(:) = 0.d0
  flap_selection(:) = 'specify'
  flap_degrees(:) = 0.d0
  weighting(:) = 1.d0

  seed_violation_handling = 'stop'
  min_thickness = 0.06d0
  max_thickness = 1000.d0
  moment_constraint_type(:) = 'use_seed'
  min_moment(:) = -1.d0
  min_te_angle = 5.d0
  check_curvature = .false.
  max_curv_reverse = 3
  curv_threshold = 0.30d0
  symmetrical = .false.
  min_flap_degrees = -5.d0
  max_flap_degrees = 15.d0

! Read operating conditions and constraints

  rewind(iunit)
  read(iunit, iostat=iostat1, nml=operating_conditions)
  if (iostat1 /= 0) then
    msg = "Unrecognized variable name in namelist "//&
    "'operating_conditions.' Check User Guide for available variables."
    call my_stop(msg)
  end if
  rewind(iunit)
  read(iunit, iostat=iostat1, nml=constraints)
  if (iostat1 /= 0) then
    msg = "Unrecognized variable name in namelist "//&
    "'constraints.' Check User Guide for available variables."
    call my_stop(msg)
  end if

! Store operating points where flap setting will be optimized

  nflap_optimize = 0
  if (use_flap .and. (.not. match_foils)) then
    do i = 1, noppoint
      if (flap_selection(i) == 'optimize') then
        nflap_optimize = nflap_optimize + 1
        flap_optimize_points(nflap_optimize) = i
      end if
    end do
  end if

! Normalize weightings for operating points

  weighting = weighting/sum(weighting(1:noppoint))

! Set default initialization options

  feasible_init = .true.
  feasible_limit = 5.0D+04
  feasible_init_attempts = 1000

! Read initialization parameters

  rewind(iunit)
  read(iunit, iostat=iostat1, nml=initialization)
  if (iostat1 /= 0) then
    msg = "Unrecognized variable name in namelist "//&
    "'initialization.' Check User Guide for available variables."
    call my_stop(msg)
  end if

! Set default particle swarm options

  pso_pop = 40
  pso_tol = 1.D-04
  pso_maxit = 300
  pso_convergence_profile = 'standard'

! Set default genetic algorithm options

  ga_pop = 80
  ga_tol = 1.D-04
  ga_maxit = 300
  parents_selection_method = 'roulette'
  parent_fraction = 0.5d0
  roulette_selection_pressure = 8.d0
  tournament_fraction = 0.1d0
  crossover_range_factor = 0.4d0
  mutant_probability = 0.3d0
  chromosome_mutation_rate = 0.01d0
  mutation_range_factor = 0.2d0

! Set default simplex search options

  simplex_tol = 1.0D-05
  simplex_maxit = 1000

  if (trim(search_type) == 'global_and_local' .or. trim(search_type) ==        &
      'global') then
  
!   Set design variables with side constraints

    if (trim(shape_functions) == 'naca') then

!     For NACA, we will only constrain the flap deflection

      allocate(constrained_dvs(nflap_optimize))
      counter = 0
      do i = nfunctions_top + nfunctions_bot + 1,                              &
             nfunctions_top + nfunctions_bot + nflap_optimize
        counter = counter + 1
        constrained_dvs(counter) = i
      end do
          
    else

!     For Hicks-Henne, also constrain bump locations and width

      allocate(constrained_dvs(2*nfunctions_top + 2*nfunctions_bot +           &
                               nflap_optimize))
      counter = 0
      do i = 1, nfunctions_top + nfunctions_bot
        counter = counter + 1
        idx = 3*(i-1) + 2      ! DV index of bump location, shape function i
        constrained_dvs(counter) = idx
        counter = counter + 1
        idx = 3*(i-1) + 3      ! Index of bump width, shape function i
        constrained_dvs(counter) = idx
      end do
      do i = 3*(nfunctions_top + nfunctions_bot) + 1,                          &
             3*(nfunctions_top + nfunctions_bot) + nflap_optimize
        counter = counter + 1
        constrained_dvs(counter) = i
      end do

    end if

    if (trim(global_search) == 'particle_swarm') then

!     Read PSO options and put them into derived type

      rewind(iunit)
      read(iunit, iostat=iostat1, nml=particle_swarm_options)
      if (iostat1 /= 0) then
        msg = "Unrecognized variable name in namelist "//&
        "'particle_swarm_options.' Check User Guide for available variables."
        call my_stop(msg)
      end if
      pso_options%pop = pso_pop
      pso_options%tol = pso_tol
      pso_options%maxspeed = initial_perturb
      pso_options%maxit = pso_maxit
      pso_options%convergence_profile = pso_convergence_profile
      pso_options%feasible_init = feasible_init
      pso_options%feasible_limit = feasible_limit
      pso_options%feasible_init_attempts = feasible_init_attempts
      pso_options%write_designs = write_designs
      if (.not. match_foils) then
        pso_options%relative_fmin_report = .true.
      else
        pso_options%relative_fmin_report = .false.
      end if

    else if (trim(global_search) == 'genetic_algorithm') then

!     Read genetic algorithm options and put them into derived type

      rewind(iunit)
      read(iunit, iostat=iostat1, nml=genetic_algorithm_options)
      if (iostat1 /= 0) then
        msg = "Unrecognized variable name in namelist "//&
        "'genetic_algorithm_options.' Check User Guide for available variables."
        call my_stop(msg)
      end if
      ga_options%pop = ga_pop
      ga_options%tol = ga_tol
      ga_options%maxit = ga_maxit
      ga_options%parents_selection_method = parents_selection_method
      ga_options%parent_fraction = parent_fraction
      ga_options%roulette_selection_pressure = roulette_selection_pressure
      ga_options%tournament_fraction = tournament_fraction
      ga_options%crossover_range_factor = crossover_range_factor
      ga_options%mutant_probability = mutant_probability
      ga_options%chromosome_mutation_rate = chromosome_mutation_rate
      ga_options%mutation_range_factor = mutation_range_factor
      ga_options%feasible_init = feasible_init
      ga_options%feasible_limit = feasible_limit
      ga_options%feasible_init_attempts = feasible_init_attempts
      ga_options%write_designs = write_designs
      if (.not. match_foils) then
        ga_options%relative_fmin_report = .true.
      else
        ga_options%relative_fmin_report = .false.
      end if

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

      rewind(iunit)
      read(iunit, iostat=iostat1, nml=simplex_options)
      if (iostat1 /= 0) then
        msg = "Unrecognized variable name in namelist "//&
        "'simplex_options.' Check User Guide for available variables."
        call my_stop(msg)
      end if
      ds_options%tol = simplex_tol
      ds_options%maxit = simplex_maxit
      ds_options%write_designs = write_designs
      if (.not. match_foils) then
        ds_options%relative_fmin_report = .true.
      else
        ds_options%relative_fmin_report = .false.
      end if

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

  rewind(iunit)
  read(iunit, iostat=iostat1, nml=xfoil_run_options)
  if (iostat1 /= 0) then
    msg = "Unrecognized variable name in namelist "//&
    "'xfoil_run_options.' Check User Guide for available variables."
    call my_stop(msg)
  end if
  rewind(iunit)
  read(iunit, iostat=iostat1, nml=xfoil_paneling_options)
  if (iostat1 /= 0) then
    msg = "Unrecognized variable name in namelist "//&
    "'xfoil_paneling_options.' Check User Guide for available variables."
    call my_stop(msg)
  end if

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
  matchfoil_file = 'none'
  read(iunit, iostat=iostat1, nml=matchfoil_options)
  if (iostat1 /= 0) then
    msg = "Unrecognized variable name in namelist "//&
    "'matchfoil_options.' Check User Guide for available variables."
    call my_stop(msg)
  end if

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
  write(*,*) " min_bump_width = ", min_bump_width
  write(*,*) " nfunctions_top = ", nfunctions_top
  write(*,*) " nfunctions_bot = ", nfunctions_bot
  write(*,*) " initial_perturb = ", initial_perturb
  write(*,*) " restart = ", restart
  write(*,*) " restart_write_freq = ", restart_write_freq
  write(*,*) " write_designs = ", write_designs
  write(*,'(A)') " /"
  write(*,*)

! Operating conditions namelist

  write(*,'(A)') " &operating_conditions"
  write(*,*) " noppoint = ", noppoint
  write(*,*) " use_flap = ", use_flap
  write(*,*) " x_flap = ", x_flap
  write(*,*) " y_flap = ", y_flap
  write(*,*)
  do i = 1, noppoint
    write(text,*) i
    text = adjustl(text)
    write(*,*) " optimization_type("//trim(text)//") = '"//                    &
               trim(optimization_type(i))//"'"
    write(*,*) " op_mode("//trim(text)//") = '"//trim(op_mode(i))//"'"
    write(*,*) " op_point("//trim(text)//") = ", op_point(i)
    write(*,'(A,es17.8)') "  reynolds("//trim(text)//") = ", reynolds(i)
    write(*,*) " mach("//trim(text)//") = ", mach(i)
    write(*,*) " flap_selection("//trim(text)//") = '"//                       &
               trim(flap_selection(i))//"'"
    write(*,*) " flap_degrees("//trim(text)//") = ", flap_degrees(i)
    write(*,*) " weighting("//trim(text)//") = ", weighting(i)
    if (i < noppoint) write(*,*)
  end do
  write(*,'(A)') " /"
  write(*,*)

! Constraints namelist

  write(*,'(A)') " &constraints"
  write(*,*) " seed_violation_handling = "//trim(seed_violation_handling)
  write(*,*) " min_thickness = ", min_thickness
  write(*,*) " max_thickness = ", max_thickness
  do i = 1, noppoint
    write(text,*) i
    text = adjustl(text)
    write(*,*) " moment_constraint_type("//trim(text)//") = "//                &
               trim(moment_constraint_type(i))
    write(*,*) " min_moment("//trim(text)//") = ", min_moment(i)
  end do
  write(*,*) " min_te_angle = ", min_te_angle
  write(*,*) " check_curvature = ", check_curvature
  write(*,*) " max_curv_reverse = ", max_curv_reverse
  write(*,*) " curv_threshold = ", curv_threshold
  write(*,*) " symmetrical = ", symmetrical
  write(*,*) " min_flap_degrees = ", min_flap_degrees
  write(*,*) " max_flap_degrees = ", max_flap_degrees
  write(*,'(A)') " /"
  write(*,*)

! Initialization namelist

  write(*,'(A)') " &initialization"
  write(*,*) " feasible_init = ", feasible_init
  write(*,*) " feasible_limit = ", feasible_limit
  write(*,*) " feasible_init_attempts = ", feasible_init_attempts
  write(*,'(A)') " /"
  write(*,*)

! Optimizer namelists

  if (trim(search_type) == 'global_and_local' .or. trim(search_type) ==        &
      'global') then

    if (trim(global_search) == 'particle_swarm') then

!     Particle swarm namelist

      write(*,'(A)') " &particle_swarm_options"
      write(*,*) " pso_pop = ", pso_options%pop
      write(*,*) " pso_tol = ", pso_options%tol
      write(*,*) " pso_maxit = ", pso_options%maxit
      write(*,*) " pso_convergence_profile = ", pso_options%convergence_profile
      write(*,'(A)') " /"
      write(*,*)

    else if (trim(global_search) == 'genetic_algorithm') then

!     Genetic algorithm options

      write(*,'(A)') " &genetic_algorithm_options"
      write(*,*) " ga_pop = ", ga_options%pop
      write(*,*) " ga_tol = ", ga_options%tol
      write(*,*) " ga_maxit = ", ga_options%maxit
      write(*,*) " parents_selection_method = ",                               &
                 ga_options%parents_selection_method
      write(*,*) " parent_fraction = ", ga_options%parent_fraction 
      write(*,*) " roulette_selection_pressure = ",                            &
                 ga_options%roulette_selection_pressure
      write(*,*) " tournament_fraction = " , ga_options%tournament_fraction
      write(*,*) " crossover_range_factor = ", ga_options%crossover_range_factor
      write(*,*) " mutant_probability = ", ga_options%mutant_probability
      write(*,*) " chromosome_mutation_rate = ",                               &
                 ga_options%chromosome_mutation_rate
      write(*,*) " mutation_range_factor = ", ga_options%mutation_range_factor
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
  if (min_bump_width <= 0.d0)                                                  &
    call my_stop("min_bump_width must be > 0.")

! Operating points

  if (noppoint < 1) call my_stop("noppoint must be > 0.")
  if ((use_flap) .and. (x_flap <= 0.0)) call my_stop("x_flap must be > 0.")
  if ((use_flap) .and. (x_flap >= 1.0)) call my_stop("x_flap must be < 1.")

  do i = 1, noppoint
    if (trim(op_mode(i)) /= 'spec-cl' .and. trim(op_mode(i)) /= 'spec-al')     &
      call my_stop("op_mode must be 'spec-al' or 'spec-cl'.")
    if (reynolds(i) <= 0.d0) call my_stop("reynolds must be > 0.")
    if (mach(i) < 0.d0) call my_stop("mach must be >= 0.")
    if (trim(flap_selection(i)) /= 'specify' .and.                             &
        trim(flap_selection(i)) /= 'optimize')                                 &
      call my_stop("flap_selection must be 'specify' or 'optimize'.")
    if (flap_degrees(i) < -90.d0) call my_stop("flap_degrees must be > -90.")
    if (flap_degrees(i) > 90.d0) call my_stop("flap_degrees must be < 90.")
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
  do i = 1, noppoint
    if (trim(moment_constraint_type(i)) /= 'use_seed' .and.                    &
      trim(moment_constraint_type(i)) /= 'specify' .and.                       &
      trim(moment_constraint_type(i)) /= 'none')                               &
      call my_stop("moment_constraint_type must be 'use_seed', 'specify', "//  &
                 "or 'none'.")
  end do
  if (min_te_angle <= 0.d0) call my_stop("min_te_angle must be > 0.")
  if (max_curv_reverse < 1) call my_stop("max_curv_reverse must be > 0.")
  if (curv_threshold <= 0.d0) call my_stop("curv_threshold must be > 0.")
  if (symmetrical)                                                             &
    write(*,*) "Mirroring top half of seed airfoil for symmetrical constraint."
  if (min_flap_degrees >= max_flap_degrees)                                    &
    call my_stop("min_flap_degrees must be less than max_flap_degrees.")
  if (min_flap_degrees <= -90.d0)                                              &
    call my_stop("min_flap_degrees must be greater than -90.")
  if (max_flap_degrees >= 90.d0)                                               &
    call my_stop("max_flap_degrees must be less than 90.")

! Initialization options
    
  if ((feasible_limit <= 0.d0) .and. feasible_init)                            &
    call my_stop("feasible_limit must be > 0.")
  if ((feasible_init_attempts < 1) .and. feasible_init)                        &
    call my_stop("feasible_init_attempts must be > 0.")

! Optimizer options

  if (trim(search_type) == 'global' .or.                                       &
       trim(search_type) == 'global_and_local') then

    if (trim(global_search) == 'particle_swarm') then

!     Particle swarm options

      if (pso_pop < 1) call my_stop("pso_pop must be > 0.")
      if (pso_tol <= 0.d0) call my_stop("pso_tol must be > 0.")
      if (pso_maxit < 1) call my_stop("pso_maxit must be > 0.")  
      if ( (trim(pso_convergence_profile) /= "standard") .and.                 &
           (trim(pso_convergence_profile) /= "exhaustive") )                   &
        call my_stop("pso_convergence_profile must be 'standard' "//&
                     "or 'exhaustive'.")

    else if (trim(global_search) == 'genetic_algorithm') then

!     Genetic algorithm options

      if (ga_pop < 1) call my_stop("ga_pop must be > 0.")
      if (ga_tol <= 0.d0) call my_stop("ga_tol must be > 0.")
      if (ga_maxit < 1) call my_stop("ga_maxit must be > 0.")
      if ( (trim(parents_selection_method) /= "roulette") .and.                &
           (trim(parents_selection_method) /= "tournament") .and.              &
           (trim(parents_selection_method) /= "random") )                      &
        call my_stop("parents_selection_method must be 'roulette', "//&
                     "'tournament', or 'random'.")
      if ( (parent_fraction <= 0.d0) .or. (parent_fraction > 1.d0) )           &
        call my_stop("parent_fraction must be > 0 and <= 1.")
      if (roulette_selection_pressure <= 0.d0)                                 &
        call my_stop("roulette_selection_pressure must be > 0.")
      if ( (tournament_fraction <= 0.d0) .or. (tournament_fraction > 1.d0) )   &
        call my_stop("tournament_fraction must be > 0 and <= 1.")
      if (crossover_range_factor < 0.d0)                                       &
        call my_stop("crossover_range_factor must be >= 0.")
      if ( (mutant_probability < 0.d0) .or. (mutant_probability > 1.d0) )      &
        call my_stop("mutant_probability must be >= 0 and <= 1.") 
      if (chromosome_mutation_rate < 0.d0)                                     &
        call my_stop("chromosome_mutation_rate must be >= 0.")
      if (mutation_range_factor < 0.d0)                                        &
        call my_stop("mutation_range_factor must be >= 0.")

    end if

  end if

  if (trim(search_type) == 'local' .or.                                        &
       trim(search_type) == 'global_and_local') then

!   Simplex options

    if (simplex_tol <= 0.d0) call my_stop("simplex_tol must be > 0.")
    if (simplex_maxit < 1) call my_stop("simplex_maxit must be > 0.")  
  
  end if

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
                                 reynolds, mach, use_flap, x_flap, y_flap,     &
                                 flap_degrees
  use airfoil_evaluation, only : xfoil_options, xfoil_geom_options
 
  character(80), intent(out) :: airfoil_file

  logical :: viscous_mode, silent_mode, fix_unconverged, reinitialize
  integer :: bl_maxit, npan
  double precision :: ncrit, xtript, xtripb, vaccel
  double precision :: cvpar, cterat, ctrrat, xsref1, xsref2, xpref1, xpref2
  integer :: i, iunit, ioerr, iostat1
  character(30) :: text

  namelist /airfoil_to_load/ airfoil_file
  namelist /operating_conditions/ noppoint, op_mode, op_point, reynolds, mach, &
            use_flap, x_flap, y_flap, flap_degrees
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

  read(iunit, iostat=iostat1, nml=airfoil_to_load)

! Set defaults for operating conditions

  noppoint = 1
  use_flap = .false.
  x_flap = 0.75d0
  y_flap = 0.d0
  op_mode(:) = 'spec-cl'
  op_point(:) = 0.d0
  reynolds(:) = 1.0D+05
  mach(:) = 0.d0
  flap_degrees = 0.d0

! Read operating conditions and constraints

  read(iunit, iostat=iostat1, nml=operating_conditions)

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

  read(iunit, iostat=iostat1, nml=xfoil_run_options)
  read(iunit, iostat=iostat1, nml=xfoil_paneling_options)

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
  write(*,*) " use_flap = ", use_flap
  write(*,*) " x_flap = ", x_flap
  write(*,*) " y_flap = ", y_flap
  do i = 1, noppoint
    write(text,*) i
    text = adjustl(text)
    write(*,*) " op_mode("//trim(text)//") = '"//trim(op_mode(i))//"'"
    write(*,*) " op_point("//trim(text)//") = ", op_point(i)
    write(*,'(A,es17.8)') "  reynolds("//trim(text)//") = ", reynolds(i)
    write(*,*) " mach("//trim(text)//") = ", mach(i)
    write(*,*) " flap_degrees("//trim(text)//") = ", flap_degrees(i)
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
