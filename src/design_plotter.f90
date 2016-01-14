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

program design_plotter

! Reads design variables from optimization and writes them as airfoil
! geometry and airfoil polars

  use vardef
  use input_output,        only : read_inputs, read_clo
  use optimization,        only : pso_options_type, ds_options_type
  use airfoil_operations,  only : get_seed_airfoil, get_split_points,          &
                                  split_airfoil, deallocate_airfoil
  use parameterization,    only : create_shape_functions,                      &
                                  deallocate_shape_functions

  implicit none

  type(airfoil_type) :: buffer_foil
  character(80) :: search_type, global_search, local_search, seed_airfoil,     &
                   airfoil_file, matchfoil_file
  character(4) :: naca_digits
  character(80) :: input_file, output_prefix
  type(pso_options_type) :: pso_options
  type(ds_options_type) :: ds_options
  integer :: pointst, pointsb
  integer, dimension(:), allocatable :: constrained_dvs
  double precision, dimension(:), allocatable :: modest, modesb

  write(*,*)
  write(*,*) 'Design plotter: writes airfoils and polars for designs generated'
  write(*,*) 'during optimization. Copyright 2014 -- 2016 Daniel Prosser.'

! Read command line arguments

  call read_clo(input_file, output_prefix)

! Read inputs from namelist file

  call read_inputs(input_file, search_type, global_search, local_search,       &
                   seed_airfoil, airfoil_file, naca_digits, nparams_top,       &
                   nparams_bot, constrained_dvs, pso_options, ds_options,      &
                   matchfoil_file)

! Load seed airfoil into memory, including transformations and smoothing

  call get_seed_airfoil(seed_airfoil, airfoil_file, naca_digits, buffer_foil)

! Split up seed airfoil into upper and lower surfaces

  call get_split_points(buffer_foil, pointst, pointsb, symmetrical)
  allocate(xseedt(pointst))
  allocate(zseedt(pointst))
  allocate(xseedb(pointsb))
  allocate(zseedb(pointsb))
  call split_airfoil(buffer_foil, xseedt, xseedb, zseedt, zseedb, symmetrical)

! Allocate memory for shape functions and create them if NACA

  if (trim(shape_functions) == 'naca') then
    allocate(modest(nparams_top))
    allocate(modesb(nparams_bot))
  else
    allocate(modest(nparams_top*3))
    allocate(modesb(nparams_bot*3))
  end if
  modest(:) = 0.d0
  modesb(:) = 0.d0
  call create_shape_functions(xseedt, xseedb, modest, modesb, shape_functions, &
                              first_time=.true.)

! Deallocate the buffer airfoil (no longer needed)

  call deallocate_airfoil(buffer_foil)

! Create airfoils and polars from designs

  call design_visualize(search_type, global_search, local_search, output_prefix)

! Deallocate other memory

  deallocate(xseedt)
  deallocate(xseedb)
  deallocate(zseedt)
  deallocate(zseedb)
  deallocate(modest)
  deallocate(modesb)
  call deallocate_shape_functions()

end program design_plotter

!=============================================================================80
!
! Driver routine to create airfoils and polars from designs
!
!=============================================================================80
subroutine design_visualize(search_type, global_search, local_search,          &
                            output_prefix)

  use vardef, only : match_foils

  implicit none
  character(*), intent(in) :: output_prefix, search_type, global_search,       &
                              local_search

  integer :: foilunit, polarunit, designnum

  designnum = 0

! Based on optimization options, select which files to read and write

  if (trim(search_type) == 'global' .or. trim(search_type) ==                  &
      'global_and_local') then

    if (trim(global_search) == 'particle_swarm') then

      call plotter(output_prefix, 'particleswarm', 'new', designnum, foilunit, &
                   polarunit)

    end if

  end if

  if (trim(search_type) == 'local') then 

    if (trim(local_search) == 'simplex') then

      call plotter(output_prefix, 'simplex', 'new', designnum, foilunit,       &
                   polarunit)

    end if

  end if

  if (trim(search_type) == 'global_and_local') then 

    if (trim(local_search) == 'simplex') then

      call plotter(output_prefix, 'simplex', 'old', designnum, foilunit,       &
                   polarunit)

    end if

  end if

! Close the files that have been written to

  close(foilunit)
  if (.not. match_foils) close(polarunit)

end subroutine design_visualize

!=============================================================================80
!
! Subroutine to create airfoils and polars from design variables
!
!=============================================================================80
subroutine plotter(writetitle, readtitle, filestat, designnum, foilunit,       &
                   polarunit)

  use vardef
  use parameterization,   only : top_shape_function, bot_shape_function,       &
                                 create_airfoil
  use airfoil_operations, only : allocate_airfoil, deallocate_airfoil
  use xfoil_driver,       only : run_xfoil, xfoil_init, xfoil_cleanup
  use airfoil_evaluation, only : xfoil_options, xfoil_geom_options

  implicit none
  character(*), intent(in) :: writetitle, readtitle, filestat
  integer, intent(inout) :: designnum
  integer, intent(out) :: foilunit, polarunit

  integer :: iunit, ioerr, nvars, designcounter
  character(100) :: readfile, foilfile, polarfile, blankchar, text
  double precision, dimension(:), allocatable :: designvars
  double precision, dimension(size(xseedt,1)) :: zt_new
  double precision, dimension(size(xseedb,1)) :: zb_new
  integer :: nmodest, nmodesb, nptt, nptb, i, dvtbnd1, dvtbnd2, dvbbnd1,       &
             dvbbnd2, flap_idx, dvcounter, oppoint, nfuncs
  double precision, dimension(noppoint) :: lift, drag, moment, viscrms
  double precision, dimension(noppoint) :: actual_flap_degrees
  double precision :: ffact

  iunit = 12
  foilunit = 13
  polarunit = 14

! Information about shape functions

  nmodest = size(top_shape_function,1)
  nmodesb = size(bot_shape_function,1)
  nptt = size(xseedt,1)
  nptb = size(xseedb,1)

! Allocate memory for airfoil

  curr_foil%npoint = size(xseedt,1) + size(xseedb) - 1
  call allocate_airfoil(curr_foil)

! Set file names 

  readfile = trim(readtitle)//'_designs.dat'
  if (writetitle == "optfoil") then      ! Probably not specified by the user
    foilfile = 'design_coordinates.dat'
    polarfile = 'design_polars.dat'
  else
    foilfile = trim(writetitle)//'_design_coordinates.dat'
    polarfile = trim(writetitle)//'_design_polars.dat'
  end if

! Open file and read number of design variables

  open(unit=iunit, file=readfile, status='old', iostat=ioerr)
  write(*,*) 'Reading designs from the file '//trim(readfile)//' ...'
  write(*,*)
  if (ioerr /= 0) then
    write(*,*) 'Error: cannot find file '//trim(readfile)
    write(*,*)
    stop
  end if
  read(iunit,'(A20,I8)') blankchar, nvars

! Allocate design variables

  allocate(designvars(nvars))

! Allocate xfoil variables

  call xfoil_init()

! Open files to write and write headers, if necessary

  if (trim(filestat) == 'new') then

    open(unit=foilunit, file=foilfile, status='replace')
    write(foilunit,'(A)') 'title="Airfoil coordinates"'
    write(foilunit,'(A)') 'variables="x" "z"'

!   Set up seed airfoil written for comparison

    designvars(:) = 0.d0
    nfuncs = nvars - nflap_optimize
    ffact = initial_perturb/(max_flap_degrees - min_flap_degrees)
    do i = nfuncs + 1, nvars
      oppoint = flap_optimize_points(i-nfuncs)
      designvars(i) = flap_degrees(oppoint)*ffact
    end do

!   Write a message about coordinates being analyzed

    write(*,*) '  Analyzing coordinates, seed airfoil ...'

!   Set modes for top and bottom surfaces

    dvtbnd1 = 1
    if (trim(shape_functions) == 'naca') then
      dvtbnd2 = nmodest
      dvbbnd2 = nmodest + nmodesb
    else
      dvtbnd2 = nmodest*3
      dvbbnd2 = nmodest*3 + nmodesb*3
    end if
    dvbbnd1 = dvtbnd2 + 1

!   Overwrite lower DVs for symmetrical airfoils (they are not used)

    if (symmetrical) then
      dvbbnd1 = 1
      dvbbnd2 = dvtbnd2
    end if

!   Create top and bottom surfaces by perturbation of seed airfoil

    call create_airfoil(xseedt, zseedt, xseedb, zseedb,                        &
                        designvars(dvtbnd1:dvtbnd2),                           &
                        designvars(dvbbnd1:dvbbnd2), zt_new, zb_new,           &
                        shape_functions, symmetrical)

!   Format seed airfoil in a single loop in derived type

    do i = 1, nptt
      curr_foil%x(i) = xseedt(nptt-i+1)
      curr_foil%z(i) = zt_new(nptt-i+1)
    end do
    do i = 1, nptb-1
      curr_foil%x(i+nptt) = xseedb(i+1)
      curr_foil%z(i+nptt) = zb_new(i+1)
    end do

!   Write coordinates to file

    write(foilunit,'(A)') 'zone t="Seed airfoil"'
    do i = 1, nptt + nptb - 1
      write(foilunit,'(2es17.8)') curr_foil%x(i), curr_foil%z(i)
    end do

!   Get actual flap angles based on design variables

    ffact = initial_perturb/(max_flap_degrees - min_flap_degrees)
    actual_flap_degrees(1:noppoint) = flap_degrees(1:noppoint)
    dvcounter = dvbbnd2 + 1
    do i = 1, nflap_optimize
      flap_idx = flap_optimize_points(i)
      actual_flap_degrees(flap_idx) = designvars(dvcounter)/ffact
      dvcounter = dvcounter + 1
    end do

!   For aerodynamic optimization, set up file for polars

    if (.not. match_foils) then

      open(unit=polarunit, file=polarfile, status='replace')
      write(polarunit,'(A)') 'title="Airfoil polars"'
      write(polarunit,'(A)') 'variables="cl" "cd"'

!     Write a message about polars being computed

      write(*,*) '  Computing polars, seed airfoil ...'

!     Run xfoil for seed airfoil

      call run_xfoil(curr_foil, xfoil_geom_options, op_point(1:noppoint),      &
                     op_mode(1:noppoint), reynolds(1:noppoint),                &
                     mach(1:noppoint), use_flap, x_flap, y_flap,               &
                     actual_flap_degrees(1:noppoint), xfoil_options, lift,     &
                     drag, moment, viscrms)

!     Write polars to file

      write(polarunit,'(A)') 'zone t="Seed airfoil polar"'
      do i = 1, noppoint
        write(polarunit,'(2es17.8)') lift(i), drag(i)
      end do

    end if

  else

    open(unit=foilunit, file=foilfile, status='old', position='append')
    if (.not. match_foils) then
      open(unit=polarunit, file=polarfile, status='old', position='append')
    end if

  end if

! Loop over designs

  do 

!   Read design variables

    read(iunit,'(A13,I8)',end=500) blankchar, designcounter
    do i = 1, nvars
      read(iunit,*) designvars(i)
    end do

!   Write a message about coordinates being analyzed

    write(text,*) designnum + designcounter
    text = adjustl(text)
    write(*,*) '  Analyzing coordinates, design '//trim(text)//' ...'

!   Set modes for top and bottom surfaces

    dvtbnd1 = 1
    if (trim(shape_functions) == 'naca') then
      dvtbnd2 = nmodest
      dvbbnd2 = nmodest + nmodesb
    else
      dvtbnd2 = nmodest*3
      dvbbnd2 = nmodest*3 + nmodesb*3
    end if
    dvbbnd1 = dvtbnd2 + 1

!   Overwrite lower DVs for symmetrical airfoils (they are not used)

    if (symmetrical) then
      dvbbnd1 = 1
      dvbbnd2 = dvtbnd2
    end if

!   Create top and bottom surfaces by perturbation of seed airfoil

    call create_airfoil(xseedt, zseedt, xseedb, zseedb,                        &
                        designvars(dvtbnd1:dvtbnd2),                           &
                        designvars(dvbbnd1:dvbbnd2), zt_new, zb_new,           &
                        shape_functions, symmetrical)

!   Format coordinates in a single loop in derived type

    do i = 1, nptt
      curr_foil%x(i) = xseedt(nptt-i+1)
      curr_foil%z(i) = zt_new(nptt-i+1)
    end do
    do i = 1, nptb-1
      curr_foil%x(i+nptt) = xseedb(i+1)
      curr_foil%z(i+nptt) = zb_new(i+1)
    end do

!   Write coordinates to file

    write(foilunit,'(A)') 'zone t="Airfoil", SOLUTIONTIME='//trim(text)
    do i = 1, nptt + nptb - 1
      write(foilunit,'(2es17.8)') curr_foil%x(i), curr_foil%z(i)
    end do

!   Get actual flap angles based on design variables

    ffact = initial_perturb/(max_flap_degrees - min_flap_degrees)
    actual_flap_degrees(1:noppoint) = flap_degrees(1:noppoint)
    dvcounter = dvbbnd2 + 1
    do i = 1, nflap_optimize
      flap_idx = flap_optimize_points(i)
      actual_flap_degrees(flap_idx) = designvars(dvcounter)/ffact
      dvcounter = dvcounter + 1
    end do

!   Run xfoil if this was an aerodynamic optimization

    if (.not. match_foils) then

!     Write a message about polars being analyzed

      write(*,*) '  Computing polars, design '//trim(text)//' ...'

      call run_xfoil(curr_foil, xfoil_geom_options, op_point(1:noppoint),      &
                     op_mode(1:noppoint), reynolds(1:noppoint),                &
                     mach(1:noppoint), use_flap, x_flap, y_flap,               &
                     actual_flap_degrees(1:noppoint), xfoil_options, lift,     &
                     drag, moment, viscrms)

!     Write polars to file

      write(polarunit,'(A)') 'zone t="Polars", SOLUTIONTIME='//trim(text)
      do i = 1, noppoint
        write(polarunit,'(2es17.8)') lift(i), drag(i)
      end do

    end if

  end do

! Close the file storing design variables

500 close(iunit)

! Write a message about output files

  write(*,*)
  write(*,*) 'Writing airfoil coordinates to the file '//trim(foilfile)//' ...'
  write(*,*)
  if (.not. match_foils) then
    write(*,*) 'Writing airfoil polars to the file '//trim(polarfile)//' ...'
    write(*,*)
  end if

! Store number of designs

  designnum = designnum + designcounter

! Deallocate design variables and airfoil

  deallocate(designvars)
  call deallocate_airfoil(curr_foil)

! Deallocate xfoil variables

  call xfoil_cleanup()

end subroutine plotter
