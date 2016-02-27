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

module genetic_algorithm

! Module containing genetic algorithm optimization routine

  implicit none

! Options type definition for genetic algorithm

  type ga_options_type

    integer :: pop                ! genetic algorithm population size
    double precision :: tol       ! tolerance in max radius of designs before
                                  !   triggering a stop condition
    integer :: maxit              ! Max steps allowed before stopping
    character(10) :: parents_selection_method 
                                  ! method for selecting designs to reproduce:
                                  !   roulette, tournament, or random
    double precision :: parent_fraction
                                  ! fraction of population selected to 
                                  !   reproduce
    double precision :: roulette_selection_pressure
                                  ! factor to increase likelihood of best
                                  !   designs being selected by roulette
    double precision :: tournament_fraction
                                  ! fraction of population considered in
                                  !   tournament selection of each parent
    double precision :: crossover_range_factor
                                  ! fraction by which parent characteristics
                                  !   can be extrapolated during crossover
    double precision :: mutant_probability
                                  ! probability of mutation occurring in an
                                  !   offspring
    double precision :: chromosome_mutation_rate
                                  ! probability of mutation in a given
                                  !   chromosome for mutants
    double precision :: mutation_range_factor
                                  ! magnitude of change in a chromosome that
                                  !   can occur during mutation, as fraction of
                                  !   xmax - xmin
    logical :: feasible_init      ! Whether to enforce initially feasible
                                  !   designs
    double precision :: feasible_limit
                                  ! Max objective function value below which
                                  !   initial designs are considered feasible
    integer :: feasible_init_attempts
                                  ! Number of attempts to try to get a feasible
                                  !   initial design
    logical :: write_designs      ! Whether to write best design each time it
                                  !   changes
    logical :: relative_fmin_report
                                  ! If .true., reports improvement over seed
                                  !   design. Otherwise, reports fmin itself.
  end type ga_options_type

  contains

!=============================================================================80
!
! Genetic algorithm optimization routine.  Recommended as a first step to 
! determine the vicinity of the global optimum, followed by a local search to
! hone in.
!
!=============================================================================80
subroutine geneticalgorithm(xopt, fmin, step, fevals, objfunc, x0, xmin, xmax, &
                            given_f0_ref, f0_ref, constrained_dvs, ga_options, &
                            restart, restart_write_freq, designcounter,        &
                            converterfunc)

  use optimization_util, only : init_random_seed, initial_designs,             &
                                design_radius, write_design 

  double precision, dimension(:), intent(inout) :: xopt
  double precision, intent(out) :: fmin
  integer, intent(out) :: step, fevals

  interface
    double precision function objfunc(x)
      double precision, dimension(:), intent(in) :: x
    end function
  end interface

  double precision, dimension(:), intent(in) :: x0, xmin, xmax
  double precision, intent(inout) :: f0_ref
  integer, dimension(:), intent(in) :: constrained_dvs
  logical, intent(in) :: given_f0_ref, restart
  type (ga_options_type), intent(in) :: ga_options
  integer, intent(in) :: restart_write_freq
  integer, intent(out) :: designcounter

  optional :: converterfunc
  interface
    integer function converterfunc(x, designcounter)
      double precision, dimension(:), intent(in) :: x
      integer, intent(in) :: designcounter
    end function
  end interface

  integer :: nconstrained, i, j, fminloc, var, stat, restartcounter
  double precision :: mincurr, f0, radius
  double precision, dimension(ga_options%pop) :: objval
  double precision, dimension(size(xmin,1),ga_options%pop) :: dv
  logical :: use_x0, converged, signal_progress

! Get f0 (reference seed design objective function)

  if (given_f0_ref) then
    f0 = f0_ref
  else
    f0 = objfunc(x0)
    f0_ref = f0
  end if

! Initialize a random seed

  call init_random_seed()

! Set up initial designs

  if (.not. restart) then
    call initial_designs(dv, objval, fevals, objfunc, xmin, xmax, use_x0, x0,  &
                         ga_options%feasible_init, ga_options%feasible_limit,  &
                         ga_options%feasible_init_attempts)
  end if

! Set up or read other initialization data

  if (.not. restart) then

!   Global best so far

    fmin = f0
    mincurr = minval(objval,1)
    fminloc = minloc(objval,1)
    xopt = dv(:,fminloc)

!   Counters

    step = 0
    designcounter = 0

  else

!   Read restart data from file

    call ga_read_restart(step, designcounter, dv, objval, fmin, xopt)
    mincurr = minval(objval,1)

  end if

! Begin optimization

  restartcounter = 1
  converged = .false.
  write(*,*) 'Genetic algorithm optimization progress:'

end subroutine geneticalgorithm

!=============================================================================80
!
! Genetic algorithm restart write routine
!
!=============================================================================80
subroutine ga_write_restart(step, designcounter, dv, objval, fmin, xopt)

  use vardef, only : output_prefix

  integer, intent(in) :: step, designcounter
  double precision, dimension(:,:), intent(in) :: dv
  double precision, dimension(:), intent(in) :: objval, xopt
  double precision, intent(in) :: fmin

  character(100) :: restfile
  integer :: iunit

!  Status notification

  restfile = 'restart_ga_'//trim(output_prefix)
  write(*,*) '  Writing genetic algorithm restart data to file '//&
             trim(restfile)//' ...'
 
! Open output file for writing

  iunit = 13
  open(unit=iunit, file=restfile, status='replace', form='unformatted')

! Write restart data

  write(iunit) step
  write(iunit) designcounter
  write(iunit) dv
  write(iunit) objval
  write(iunit) fmin
  write(iunit) xopt

! Close restart file

  close(iunit)

! Status notification

  write(*,*) '  Successfully wrote genetic algorithm restart file.'

end subroutine ga_write_restart

!=============================================================================80
!
! Genetic algorithm restart read routine
!
!=============================================================================80
subroutine ga_read_restart(step, designcounter, dv, objval, fmin, xopt)

  use vardef, only : output_prefix

  integer, intent(inout) :: step, designcounter
  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objval, xopt
  double precision, intent(out) :: fmin

  character(100) :: restfile
  integer :: iunit, ioerr

! Status notification

  restfile = 'restart_ga_'//trim(output_prefix)
  write(*,*) '  Reading genetic algorithm restart data from file '//&
             trim(restfile)//' ...'
 
! Open output file for reading

  iunit = 13
  open(unit=iunit, file=restfile, status='old', form='unformatted',            &
       iostat=ioerr)
  if (ioerr /= 0) then
    write(*,*) 'Error: could not find input file '//trim(restfile)//'.'
    write(*,*)
    stop
  end if

! Read restart data

  read(iunit) step
  read(iunit) designcounter
  read(iunit) dv
  read(iunit) objval
  read(iunit) fmin
  read(iunit) xopt

! Close restart file

  close(iunit)

! Status notification

  write(*,*) 'Successfully read genetic algorithm restart file.'
  write(*,*)

end subroutine ga_read_restart

!=============================================================================80
!
! Selects unique parent designs for reproduction
!
!=============================================================================80
subroutine parents_selection(objval, method, beta, tournament_fraction,        &
                             idxparents)

  use optimization_util, only : pop_double_vector, pop_integer_vector

  double precision, dimension(:), intent(in) :: objval
  character(*), intent(in) :: method
  double precision, intent(in) :: beta, tournament_fraction
  integer, dimension(:), intent(inout) :: idxparents

  integer :: i, ndesigns, nparents, nconsidered, idx
  integer, dimension(size(objval,1)) :: idxconsidered
  double precision, dimension(size(objval,1)) :: objvalconsidered

! Initialize vectors

  ndesigns = size(objval,1)
  nparents = size(idxparents,1)
  objvalconsidered = objval
  do i = 1, ndesigns
    idxconsidered(i) = i
  end do

! Select nparents designs

  nconsidered = ndesigns

  do i = 1, nparents

!   Select a single design as a parent

    call single_parent_selection(objvalconsidered, nconsidered, method, beta, &
                                 tournament_fraction, idx)
    idxparents(i) = idxconsidered(idx)

!   Pop selected parent out of vectors

    call pop_integer_vector(idxconsidered, nconsidered, idx)
    call pop_double_vector(objvalconsidered, nconsidered, idx)
    nconsidered = nconsidered - 1

  end do

end subroutine parents_selection

!=============================================================================80
!
! Selects a single design for reproduction
!
!=============================================================================80
subroutine single_parent_selection(objvalconsidered, nconsidered, method, beta,&
                                   tournament_fraction, idx)

  double precision, dimension(:), intent(in) :: objvalconsidered
  integer, intent(in) :: nconsidered
  character(*), intent(in) :: method
  double precision, intent(in) :: beta, tournament_fraction
  integer, intent(out) :: idx

  continue

end subroutine single_parent_selection

end module genetic_algorithm
