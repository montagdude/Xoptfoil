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

module optimization

! Module containing optimization routines

  implicit none

! Options type definition for PSO

  type pso_options_type

    integer :: pop                ! particle swarm population size
    double precision :: speed_tol ! tolerance in particle speed, as fraction of
                                  !   initial max speed, before triggering a
                                  !   stop condition
    double precision :: maxspeed  ! Max speed allowed for particles
    integer :: maxit              ! Max steps allowed before stopping
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
    character(10) :: convergence_profile
                                  ! 'standard' or 'exhaustive'; exhaustive takes
                                  !   longer but finds better solutions 
  end type pso_options_type

! Options type for direct searches

  type ds_options_type

    double precision :: tol       ! tolerance in simplex diameter before
                                  !   triggering a stop condition
    integer :: maxit              ! Max steps allowed before stopping
    logical :: write_designs      ! Whether to write best design each time it
                                  !   changes
    logical :: relative_fmin_report
                                  ! If .true., reports improvement over seed
                                  !   design. Otherwise, reports fmin itself.
  end type ds_options_type

  contains

!=============================================================================80
!
! Particle swarm optimization routine. Recommended as a first step to determine
! the vicinity of the global optimum, followed by a gradient-based search to
! refine the optimization.
!
!=============================================================================80
subroutine particleswarm(xopt, fmin, step, fevals, objfunc, x0, xmin, xmax,    &
                         given_f0_ref, f0_ref, constrained_dvs, pso_options,   &
                         converterfunc)

  use math_deps,          only : norm_2

! The following are only needed if doing xfoil airfoil optimization

#ifdef airfoil_optimization

  use xfoil_driver,       only : xfoil_init, xfoil_cleanup
  use vardef,             only : nparams_top, nparams_bot, shape_functions,    &
                                 xseedt, xseedb, curr_foil
  use parameterization,   only : create_shape_functions,                       &
                                 deallocate_shape_functions
  use airfoil_operations, only : allocate_airfoil, deallocate_airfoil

  double precision, dimension(:), allocatable :: modest, modesb

#endif

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
  logical, intent(in) :: given_f0_ref
  type (pso_options_type), intent(in) :: pso_options

  optional :: converterfunc
  interface
    double precision function converterfunc(x)
      double precision, dimension(:), intent(in) :: x
    end function
  end interface

  integer :: nvars, nconstrained, i, j, fminloc, designcounter, var
  double precision :: c1, c2, whigh, wlow, convrate, maxspeed, wcurr, mincurr, &
                      f0, stat
  double precision, dimension(:), allocatable :: objval, minvals, randvec1,    &
                                                 randvec2, speed
  double precision, dimension(:,:), allocatable :: dv, vel, bestdesigns
  logical :: use_x0, converged

  nvars = size(xmin,1)
  nconstrained = size(constrained_dvs,1)

! PSO tuning variables

  if (trim(pso_options%convergence_profile) == "standard") then

    c1 = 1.2d0         ! particle-best trust factor
    c2 = 1.2d0         ! swarm-best trust factor
    whigh = 1.4d0      ! starting inertial parameter
    wlow = 0.6d0       ! ending inertial parameter
    convrate = 0.05d0  ! inertial parameter reduction rate

  else if (trim(pso_options%convergence_profile) == "exhaustive") then

    c1 = 1.4d0         ! particle-best trust factor
    c2 = 1.0d0         ! swarm-best trust factor
    whigh = 1.8d0      ! starting inertial parameter
    wlow = 0.8d0       ! ending inertial parameter
    convrate = 0.02d0  ! inertial parameter reduction rate

  else
    write(*,*) "Error in particleswarm: convergence mode should be"//          &
               "'standard' or 'exhaustive'."
    stop
  end if

! Speed limits

  maxspeed = abs(pso_options%maxspeed)
  if (maxspeed > maxval(xmax - xmin)) then
    maxspeed = maxval(xmax - xmin)
  elseif (maxspeed < 1.0D-14) then
    maxspeed = maxval(xmax - xmin)
  end if

! Memory allocation and initialization

  allocate(dv(nvars,pso_options%pop))
  allocate(vel(nvars,pso_options%pop))
  allocate(objval(pso_options%pop))
  allocate(bestdesigns(nvars,pso_options%pop))
  allocate(minvals(pso_options%pop))
  allocate(speed(pso_options%pop))
  allocate(randvec1(nvars))
  allocate(randvec2(nvars))
  dv(:,:) = 0.d0
  vel(:,:) = 0.d0
  objval(:) = 0.d0

!$omp parallel default(shared) private(i, j)

! To allocate private memory for airfoil optimization on each thread

#ifdef airfoil_optimization

! Allocate memory for shape functions

!$omp master
  if (trim(shape_functions) == 'naca') then
    allocate(modest(nparams_top))
    allocate(modesb(nparams_bot))
  else
    allocate(modest(nparams_top*3))
    allocate(modesb(nparams_bot*3))
  end if
  modest(:) = 0.d0
  modesb(:) = 0.d0
!$omp end master
!$omp barrier

! For NACA, this will create the shape functions.  For Hicks-Henne,
! it will just allocate them.

  call create_shape_functions(xseedt, xseedb, modest, modesb,                  &
                              shape_functions, first_time=.true.)

! Allocate memory for working airfoil on each thread

  curr_foil%npoint = size(xseedt,1) + size(xseedb,1) - 1
  call allocate_airfoil(curr_foil)

! Allocate memory for xfoil

  call xfoil_init()

#endif

! Get f0 (reference seed design objective function)

!$omp master
  if (given_f0_ref) then
    f0 = f0_ref
  else 
    f0 = objfunc(x0)
    f0_ref = f0
  end if
!$omp end master
!$omp barrier

! Initialize a random seed

  call init_random_seed()

! Set up initial designs

  use_x0 = .true.
  call initial_designs(dv, objval, fevals, objfunc, xmin, xmax, use_x0, x0,    &
                       pso_options%feasible_init, pso_options%feasible_limit,  &
                       pso_options%feasible_init_attempts)

!$omp master

! Initial velocities which may be positive or negative

  call random_number(vel)
  vel = 2.d0*maxspeed*(vel - 0.5d0)

! Matrix of best designs for each particle and vector of their values

  bestdesigns = dv
  minvals = objval

! Global best so far

  fmin = minval(objval,1)
  fminloc = minloc(objval,1)
  xopt = dv(:,fminloc)

! Write design to file if requested
! converterfunc is an optional function supplied to convert design variables 
!   into something more useful.  If not supplied, the design variables
!   themselves are written to a file.

  if (pso_options%write_designs) then
    designcounter = 1
    if (present(converterfunc)) then
      stat = converterfunc(xopt)
    else
      call write_design('particleswarm_designs.dat', 'new', xopt, designcounter)
    end if
    designcounter = designcounter + 1
  end if

! Begin optimization

  converged = .false.
  step = 0
  wcurr = whigh
  write(*,*) 'Particle swarm optimization progress:'

!$omp end master
!$omp barrier

  optimization_loop: do while (.not. converged)

!$omp master

!   Increase iteration counter

    step = step + 1

!$omp end master
!$omp barrier

!$omp do

!   Update each particle's position, evaluate objective function, etc.

    do i = 1, pso_options%pop

!     Impose speed limit

      if (speed(i) > maxspeed) then
        vel(:,i) = maxspeed*vel(:,i)/speed(i)
      end if

!     Update position and bring back to side constraints if necessary

      dv(:,i) = dv(:,i) + vel(:,i)
      do j = 1, nconstrained
        var = constrained_dvs(j)
        if (dv(var,i) < xmin(var)) then
          dv(var,i) = xmin(var)
          call random_number(speed(i))
          vel(var,i) = -speed(i)*vel(var,i)
        elseif (dv(var,i) > xmax(var)) then
          dv(var,i) = xmax(var)
          call random_number(speed(i))
          vel(var,i) = -speed(i)*vel(var,i)
        end if
      end do

!     Evaluate objective function and update local best design if appropriate

      objval(i) = objfunc(dv(:,i))
      if (objval(i) < minvals(i)) then
        minvals(i) = objval(i)
        bestdesigns(:,i) = dv(:,i)
      end if

    end do

!$omp end do

!$omp master

!   Update best overall design, if appropriate

    mincurr = minval(objval,1)
    fminloc = minloc(objval,1)
    if (mincurr < fmin) then

      xopt = dv(:,fminloc)
      fmin = mincurr

!     Write design to file if requested
!     converterfunc is an optional function supplied to convert design variables
!       into something more useful.  If not supplied, the design variables
!       themselves are written to a file.

      if (pso_options%write_designs) then
        if (present(converterfunc)) then
          stat = converterfunc(xopt)
        else
          call write_design('particleswarm_designs.dat', 'old', xopt,          &
                            designcounter)
        end if
        designcounter = designcounter + 1
      end if
    end if

!$omp end master
!$omp barrier

!$omp do

!   Update velocity of each particle

    do i = 1, pso_options%pop
      call random_number(randvec1)
      call random_number(randvec2)
      vel(:,i) = wcurr*vel(:,i) + c1*randvec1*(bestdesigns(:,i) - dv(:,i)) +   &
                                  c2*randvec2*(xopt - dv(:,i))
      speed(i) = norm2(vel(:,i))
    end do

!$omp end do

!$omp master

!   Reduce inertial parameter

    wcurr = wcurr - convrate*(wcurr - wlow)

!   Display progress

    if (pso_options%relative_fmin_report) then
      write(*,*) '  Iteration: ', step, '  % Improvement over seed: ',         &
                 (f0 - fmin)/f0*100
    else
      write(*,*) '  Iteration: ', step, ' Minimum objective function value: ', &
                 fmin
    end if
    
!   Evaluate convergence

    if ( (maxval(speed) > pso_options%speed_tol*maxspeed) .and.                &
         (step < pso_options%maxit) ) then
      converged = .false.
    else
      converged = .true.
      if (step == pso_options%maxit) then
        write(*,*) 'Warning: PSO optimizer forced to exit due to the max number'
        write(*,*) '         of iterations being reached.'
      end if
    end if 

!$omp end master
!$omp barrier

  end do optimization_loop

! To deallocate memory for airfoil optimization on each thread

#ifdef airfoil_optimization

!$omp master
  deallocate(modest)
  deallocate(modesb)
!$omp end master
  call deallocate_shape_functions()
  call deallocate_airfoil(curr_foil)
  call xfoil_cleanup()

#endif

!$omp end parallel

! Calculate number of function evaluations
      
  fevals = fevals + step*pso_options%pop

! Memory deallocation

  deallocate(dv)
  deallocate(vel)
  deallocate(objval)
  deallocate(bestdesigns)
  deallocate(minvals)
  deallocate(speed)
  deallocate(randvec1)
  deallocate(randvec2)

end subroutine particleswarm

!=============================================================================80
!
! Nelder-Mead simplex search algorithm
!
!=============================================================================80
subroutine simplex_search(xopt, fmin, step, fevals, objfunc, x0, given_f0_ref, &
                          f0_ref, ds_options, converterfunc)

! The following are only needed if doing xfoil airfoil optimization

#ifdef airfoil_optimization

  use xfoil_driver,       only : xfoil_init, xfoil_cleanup
  use vardef,             only : nparams_top, nparams_bot, shape_functions,    &
                                 xseedt, xseedb, curr_foil 
  use parameterization,   only : create_shape_functions,                       &
                                 deallocate_shape_functions
  use airfoil_operations, only : allocate_airfoil, deallocate_airfoil

  double precision, dimension(:), allocatable :: modest, modesb

#endif

  double precision, dimension(:), intent(inout) :: xopt
  double precision, intent(out) :: fmin
  integer, intent(out) :: step, fevals

  interface
    double precision function objfunc(x)
      double precision, dimension(:), intent(in) :: x
    end function
  end interface

  double precision, dimension(:), intent(in) :: x0
  double precision, intent(inout) :: f0_ref
  logical, intent(in) :: given_f0_ref
  type (ds_options_type), intent(in) :: ds_options

  optional :: converterfunc
  interface
    double precision function converterfunc(x)
      double precision, dimension(:), intent(in) :: x
    end function
  end interface

  double precision, dimension(size(x0,1),size(x0,1)+1) :: dv
  double precision, dimension(size(x0,1)+1) :: objvals
  double precision, dimension(size(x0,1)) :: xcen, xr, xe, xc

  double precision :: rho, xi, gam, sigma, fr, fe, fc, dist, diam, errval, f0, &
                      stat
  integer :: i, j, k, nvars, designcounter, nsame
  logical :: converged, needshrink
  character(3) :: filestat

! Standard Nelder-Mead constants

  rho = 1.d0
  xi = 2.d0
  gam = 0.5d0
  sigma = 0.5d0

! To allocate private memory for airfoil optimization on each thread

#ifdef airfoil_optimization

! Allocate memory for shape functions

  if (trim(shape_functions) == 'naca') then
    allocate(modest(nparams_top))
    allocate(modesb(nparams_bot))
  else
    allocate(modest(nparams_top*3))
    allocate(modesb(nparams_bot*3))
  end if
  modest(:) = 0.d0
  modesb(:) = 0.d0

! For NACA, this will create the shape functions.  For Hicks-Henne,
! it will just allocate them.

  call create_shape_functions(xseedt, xseedb, modest, modesb,                  &
                              shape_functions, first_time=.true.)

! Allocate memory for working airfoil on each thread

  curr_foil%npoint = size(xseedt,1) + size(xseedb,1) - 1
  call allocate_airfoil(curr_foil)

! Allocate memory for xfoil

  call xfoil_init()

#endif

! Get f0 (reference seed design objective function)

  if (given_f0_ref) then
    f0 = f0_ref
  else 
    f0 = objfunc(x0)
    f0_ref = f0
  end if

! Set up initial simplex

  fevals = 0
  nvars = size(x0,1)
  do j = 1, nvars
    do i = 1, nvars
      if (i == j) then
        if (x0(i) == 0.d0) then
          dv(i,j) = 0.00025d0
        else
          dv(i,j) = 1.05d0*x0(i)
        end if
      else
        dv(i,j) = x0(i)
      end if
    end do
    objvals(j) = objfunc(dv(:,j))
    fevals = fevals + 1
  end do

  dv(:,nvars+1) = x0
  objvals(nvars+1) = objfunc(x0)
  fevals = fevals + 1
  fmin = minval(objvals)

! Iterative procedure for optimization
 
  step = 0
  needshrink = .false.
  converged = .false.
  designcounter = 1
  nsame = 0
  write(*,*) 'Simplex optimization progress:'
  main_loop: do while (.not. converged)

    step = step + 1
    if (step == ds_options%maxit) converged = .true.
    
!   Sort according to ascending objective function value

    call bubble_sort(dv, objvals)
    errval = abs((objvals(1) - fmin)/fmin)
    if (errval < 1.D-15) then
      nsame = nsame + 1
    else
      nsame = 0
    end if
    fmin = objvals(1)

!   Write design to file if requested
!   converterfunc is an optional function supplied to convert design variables
!     into something more useful.  If not supplied, the design variables
!     themselves are written to a file.

    if (ds_options%write_designs .and. designcounter == 1) then
      filestat = 'new'
    else 
      filestat = 'old'
    end if

    if (ds_options%write_designs) then
      if (present(converterfunc)) then
        stat = converterfunc(dv(:,1))
      else
        call write_design('simplex_designs.dat', filestat, dv(:,1),            &
                          designcounter)
      end if
      designcounter = designcounter + 1
    end if

!   Compute max distance between simplex vertices (designs)

    diam = 0.d0
    do j = 2, nvars + 1
      dist = 0.d0
      do k = 1, nvars
        dist = dist + (dv(k,1) - dv(k,j))**2.d0
      end do
      dist = sqrt(dist)
      if (dist > diam) diam = dist
    end do

!   Check for convergence

    if (diam < ds_options%tol) converged = .true.

!   Display progress

    if (ds_options%relative_fmin_report) then
      write(*,*) '  Iteration: ', step, '  % Improvement over seed: ',         &
                 (f0 - fmin)/f0*100
    else
      write(*,*) '  Iteration: ', step, ' Minimum objective function value: ', &
                 fmin
    end if

!   Compute the centroid of the best nvals designs

    xcen(:) = 0.d0
    do i = 1, nvars
      xcen = xcen + dv(:,i)
    end do
    xcen = xcen/dble(nvars)

!   Compute the reflection point and evaluate its objective function value

    xr = (1.d0 + rho)*xcen - rho*dv(:,nvars+1)
    fr = objfunc(xr)
    fevals = fevals + 1

    expand_or_contract: if (objvals(1) <= fr .and. fr < objvals(nvars)) then

!      Accept reflection point

       dv(:,nvars+1) = xr
       objvals(nvars+1) = fr
       cycle

    elseif (fr < objvals(1)) then

!     Expand

      xe = (1.d0 + rho*xi)*xcen - rho*xi*dv(:,nvars+1)
      fe = objfunc(xe)
      fevals = fevals + 1
      if (fe < fr) then
        dv(:,nvars+1) = xe
        objvals(nvars+1) = fe
      else
        dv(:,nvars+1) = xr
        objvals(nvars+1) = fr
      end if
      cycle

    elseif (fr >= objvals(nvars)) then

!     Outside contraction

        contraction: if (fr < objvals(nvars+1)) then

          xc = (1.d0 + rho*gam)*xcen - rho*gam*dv(:,nvars+1)
          fc = objfunc(xc)
          fevals = fevals + 1

          if (fc < fr) then
            dv(:,nvars+1) = xc
            objvals(nvars+1) = fc
            needshrink = .false.
          else
            needshrink = .true.
          end if

!       Inside contraction

        else 

          xc = (1.d0 - gam)*xcen + gam*dv(:,nvars+1)
          fc = objfunc(xc)
          fevals = fevals + 1
          
          if (fc < objvals(nvars+1) ) then
            dv(:,nvars+1) = xc
            objvals(nvars+1) = fc
            needshrink = .false.
          else
            needshrink = .true.
          end if

        end if contraction

!       Shrink

        shrink: if (needshrink) then

          do i = 2, nvars + 1
            dv(:,i) = dv(:,1) + sigma*(dv(:,i) - dv(:,1))
            objvals(i) = objfunc(dv(:,i))
            fevals = fevals + 1
          end do
          cycle

        else

          cycle

        end if shrink

    end if expand_or_contract

  end do main_loop

! Sort one more time according to ascending objective function value

  call bubble_sort(dv, objvals)
  xopt = dv(:,1)
  fmin = objvals(1)

! Compute max distance between simplex vertices (designs)

  diam = 0.d0
  do j = 2, nvars + 1
    dist = 0.d0
    do k = 1, nvars
      dist = dist + (dv(k,1) - dv(k,j))**2.d0
    end do
    dist = sqrt(dist)
    if (dist > diam) diam = dist
  end do

! Display warning if max iterations are reached
  
  if (step == ds_options%maxit .and. (diam >= ds_options%tol)) then
    write(*,*) 'Warning: Simplex optimizer forced to exit due to the max number'
    write(*,*) '         of iterations being reached.'
  end if

#ifdef airfoil_optimization

  deallocate(modest)
  deallocate(modesb)
  call deallocate_shape_functions()
  call deallocate_airfoil(curr_foil)
  call xfoil_cleanup()

#endif

end subroutine simplex_search

!=============================================================================80
!
! Initializes a random seed (subroutine from gcc.gnu.org)
!
!=============================================================================80
subroutine init_random_seed()

! For ifort compatibility
#ifdef intel_compilers
  use ifport, only : getpid  
#endif

  integer, dimension(:), allocatable :: myseed
  integer :: i, n, un, istat, dt(8), pid, t(2), s
  integer(8) :: count, tms
  
  call random_seed(size = n)
  allocate(myseed(n))

  ! First try if the OS provides a random number generator

  open(newunit=un, file="/dev/urandom", access="stream",                       &
       form="unformatted", action="read", status="old", iostat=istat)

  if (istat == 0) then

     read(un) myseed
     close(un)

  else

     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.

     call system_clock(count)
     if (count /= 0) then
        t = transfer(count, t)
     else
        call date_and_time(values=dt)
        tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24 * 60 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
        t = transfer(tms, t)
     end if

     s = ieor(t(1), t(2))
     pid = getpid() + 1099279 ! Add a prime
     s = ieor(s, pid)
     if (n >= 3) then
        myseed(1) = t(1) + 36269
        myseed(2) = t(2) + 72551
        myseed(3) = pid
        if (n > 3) then
           myseed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
        end if
     else
        myseed = s + 37 * (/ (i, i = 0, n - 1 ) /)
     end if

  end if

  call random_seed(put=myseed)
  deallocate(myseed)

end subroutine init_random_seed

!=============================================================================80
!
! Creates initial designs and tries to make them feasible, if desired
!
!=============================================================================80
subroutine initial_designs(dv, objval, fevals, objfunc, xmin, xmax, use_x0,    &
                           x0, feasible_init, feasible_limit, attempts)

  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objval
  integer, intent(out) :: fevals
  double precision, dimension(:), intent(in) :: xmin, xmax, x0
  logical, intent(in) :: use_x0, feasible_init
  double precision, intent(in) :: feasible_limit
  integer, intent(in) :: attempts

  interface
    double precision function objfunc(x)
      double precision, dimension(:), intent(in) :: x
    end function
  end interface

  integer :: i, pop, nvars, initcount
  character(30) :: text1, text2, text3
  double precision, dimension(:), allocatable :: randvec1, designstore
  double precision :: minstore

! Initial settings and memory allocation

  nvars = size(dv,1)
  pop = size(dv,2)
  allocate(randvec1(nvars))
  allocate(designstore(nvars))

!$omp master

  fevals = pop

! Set up matrix of random numbers for initial designs
  
  call random_number(dv)

! Initial population of designs set between xmin and xmax

  write(*,*) 'Generating and evaluating initial designs ...'
  write(*,*)

!$omp end master
!$omp barrier

  if (use_x0) then
    dv(:,1) = x0
    objval(1) = objfunc(x0)
!$omp do
    do i = 2, pop
      dv(:,i) = maxval(xmax - xmin)*dv(:,i) + xmin
      objval(i) = objfunc(dv(:,i))
    end do
!$omp end do
  else
!$omp do
    do i = 1, pop
      dv(:,i) = maxval(xmax - xmin)*dv(:,i) + xmin
      objval(i) = objfunc(dv(:,i))
    end do
!$omp end do
  end if

! Enforce initially feasible designs

  if (feasible_init) then

    write(text1,*) attempts
    text1 = adjustl(text1)
!$omp master
    write(*,*) 'Checking initial designs for feasibility ...'
    write(*,*) '  (using a max of '//trim(text1)//' initialization attempts)'
    write(*,*)
!$omp end master


!$omp do
    do i = 1, pop

      write(text2,*) i
      text2 = adjustl(text2)
      initcount = 0
      minstore = objval(i)
      designstore = dv(:,i)

!     Take a number of tries to fix infeasible designs

      do while ((initcount <= attempts) .and.                                  &
               (objval(i) >= feasible_limit))
        call random_number(randvec1)
        dv(:,i) = maxval(xmax - xmin)*randvec1 + xmin
        objval(i) = objfunc(dv(:,i))
        if (objval(i) < minstore) then
          minstore = objval(i)
          designstore = dv(:,i)
        end if
        initcount = initcount + 1
!$omp critical
        fevals = fevals + 1
!$omp end critical
      end do

!     Pick the best design tested if a feasible design was not found

      if ((initcount > attempts) .and. (objval(i) >= feasible_limit)) then
        dv(:,i) = designstore
        objval(i) = minstore
      end if

!     Write a message about the feasibility of initial designs

      if ((initcount > attempts) .and. (objval(i) >= feasible_limit)) then
        write(*,*) ' Design '//trim(text2)//' is infeasible and was not'//     &
                   ' fixed within '//trim(text1)//' reinitialization attempts.'
      elseif ((initcount <= attempts) .and. (initcount > 0) .and.              &
              (objval(i) < feasible_limit)) then
        write(text3,*) initcount
        text3 = adjustl(text3)
        write(*,*) ' Design '//trim(text2)//' was initially infeasible but'//  &
                   ' was fixed after '//trim(text3)//                          &
                   ' reinitialization attempts.'
      else
        write(*,*) ' Design '//trim(text2)//' is feasible.' 
      end if

    end do

!$omp end do

!$omp master
    write(*,*)
!$omp end master

  end if

! Memory deallocation

  deallocate(randvec1)
  deallocate(designstore)

end subroutine initial_designs

!=============================================================================80
!
! Sorts a set of designs according to their objective function value
!
!=============================================================================80
subroutine bubble_sort(dv, objvals)

  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objvals

  double precision, dimension(size(dv,1),size(dv,2)) :: tempdv
  double precision, dimension(size(dv,2)) :: tempvals
  integer, dimension(size(dv,2)) :: finalorder, temporder
  integer :: nvars, ndesigns, i, sortcounter
  logical :: sorted

  nvars = size(dv,1)
  ndesigns = size(dv,2)

! Set up indexing array

  do i = 1, ndesigns
    finalorder(i) = i
  end do
  temporder = finalorder

! Bubble sorting algorithm

  sorted = .false.
  tempvals = objvals
  do while (.not. sorted)

    sortcounter = 0
    do i = 1, ndesigns - 1
      if (objvals(i+1) < objvals(i)) then

!       Flip the order of these elements. temp arrays are to preserve values.

        tempvals(i) = objvals(i+1)
        tempvals(i+1) = objvals(i)
        temporder(i) = finalorder(i+1)
        temporder(i+1) = finalorder(i)
        finalorder(i:i+1) = temporder(i:i+1)
        objvals(i:i+1) = tempvals(i:i+1)
        sortcounter = sortcounter + 1

      end if
    end do
    if (sortcounter == 0) sorted = .true.
    
  end do

! Use indexing array to rearrange order of designs

  do i = 1, ndesigns
    tempdv(:,i) = dv(:,finalorder(i))
  end do
  dv = tempdv

end subroutine bubble_sort

!=============================================================================80
!
! Writes design variables to file
!
!=============================================================================80
subroutine write_design(filename, filestat, variables, counter)

  character(*), intent(in) :: filename, filestat
  double precision, dimension(:), intent(in) :: variables
  integer, intent(in) :: counter

  integer, save :: iunit
  integer :: nvars, i
  character(30) :: text

  nvars = size(variables,1)
  iunit = 17

! Open the file and write to it if requested

  if (trim(filestat) == 'new') then
    open(unit=iunit, file=filename, status='replace')
    write(text,*) nvars
    text = adjustl(text)
    write(iunit,'(A)') 'Number of variables: '//trim(text)
  else
    open(unit=iunit, file=filename, status='old', position='append')
  end if

! Write iteration number and the design variables to file

  write(text,*) counter
  text = adjustl(text)
  write(iunit,'(A)') 'Design number '//trim(text)
  do i = 1, nvars
    write(iunit,'(es25.16)') variables(i)
  end do

! Close the file 

  close(iunit)

end subroutine write_design

end module optimization
