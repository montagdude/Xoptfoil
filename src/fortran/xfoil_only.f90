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

program xfoil_only

! Runs xfoil by for an airfoil, but still uses the same input file as xoptfoil
! Doesn't do any transformations or scaling on the input airfoil

  use vardef
  use input_output,       only : read_inputs, read_clo, choose_airfoil
  use particle_swarm,     only : pso_options_type
  use genetic_algorithm,  only : ga_options_type
  use simplex_search,     only : ds_options_type
  use airfoil_evaluation, only : xfoil_options, xfoil_geom_options
  use airfoil_operations, only : load_airfoil, naca_four_digit,                &
                                 deallocate_airfoil
  use xfoil_driver,       only : run_xfoil, xfoil_init, xfoil_cleanup

  implicit none

  type(airfoil_type) :: foil
  character(80) :: search_type, global_search, local_search, seed_airfoil,     &
                   airfoil_file, matchfoil_file
  character(4) :: naca_digits
  character(80) :: input_file
  type(pso_options_type) :: pso_options
  type(ga_options_type) :: ga_options
  type(ds_options_type) :: ds_options
  integer, dimension(:), allocatable :: constrained_dvs
  integer :: restart_write_freq
  logical :: restart
  double precision, dimension(:), allocatable :: alpha, lift, drag, moment,    &
                                                 viscrms

! Set default names and read command line arguments

  input_file = 'inputs.txt'
  output_prefix = 'optfoil'
  call read_clo(input_file, output_prefix)

! Read inputs from namelist file

  call read_inputs(input_file, search_type, global_search, local_search,       &
                   seed_airfoil, airfoil_file, naca_digits, nparams_top,       &
                   nparams_bot, restart, restart_write_freq, constrained_dvs,  &
                   pso_options, ga_options, ds_options, matchfoil_file)
  xfoil_options%silent_mode = .false. 

! Allocate some things

  allocate(alpha(noppoint))
  allocate(lift(noppoint))
  allocate(drag(noppoint))
  allocate(moment(noppoint))
  allocate(viscrms(noppoint))

! Ask which airfoil to analyze

  call choose_airfoil(seed_airfoil, airfoil_file, naca_digits)

! Get airfoil to analyze, but don't do any transformations

  if (trim(seed_airfoil) == "from_file") then
    call load_airfoil(airfoil_file, foil)
  else if (trim(seed_airfoil) == "four_digit") then
    call naca_four_digit(naca_digits, 200, foil)
  end if

! Allocate xfoil variables

  call xfoil_init()

! Run xfoil

  call run_xfoil(foil, xfoil_geom_options, op_point(1:noppoint),               &
                 op_mode(1:noppoint), reynolds(1:noppoint), mach(1:noppoint),  &
                 use_flap, x_flap, y_flap, flap_degrees(1:noppoint),           &
                 xfoil_options, alpha, lift, drag, moment, viscrms)

! Deallocate xfoil variables

  call xfoil_cleanup()

! Deallocate some things

  call deallocate_airfoil(foil)
  deallocate(alpha)
  deallocate(lift)
  deallocate(drag)
  deallocate(moment)
  deallocate(viscrms)

end program xfoil_only
