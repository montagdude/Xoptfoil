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

module optimization_util

! Module containing optimization routines

  implicit none

  contains

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
        finalorder(i) = temporder(i)
        finalorder(i+1) = temporder(i+1)
        objvals(i) = tempvals(i)
        objvals(i+1) = tempvals(i+1)
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

end module optimization_util
