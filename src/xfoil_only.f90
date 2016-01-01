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

  use vardef,             only : airfoil_type, max_op_points, noppoint,        &
                                 op_mode, op_point, reynolds, mach, use_flap,  &
                                 x_flap, y_flap, flap_degrees
  use input_output,       only : read_inputs_xfoil_only
  use airfoil_evaluation, only : xfoil_options, xfoil_geom_options
  use airfoil_operations, only : load_airfoil, deallocate_airfoil
  use xfoil_driver,       only : run_xfoil, xfoil_init, xfoil_cleanup

  type(airfoil_type) :: foil
  character(80) :: airfoil_file
  double precision, dimension(:), allocatable :: lift, drag, moment, viscrms

! Read inputs from namelist file

  call read_inputs_xfoil_only(airfoil_file)

! Allocate some things

  allocate(lift(noppoint))
  allocate(drag(noppoint))
  allocate(moment(noppoint))
  allocate(viscrms(noppoint))

! Load airfoil from file

  call load_airfoil(airfoil_file, foil)

! Allocate xfoil variables

  call xfoil_init()

! Run xfoil

  call run_xfoil(foil, xfoil_geom_options, op_point(1:noppoint),               &
                 op_mode(1:noppoint), reynolds(1:noppoint), mach(1:noppoint),  &
                 use_flap, x_flap, y_flap, flap_degrees(1:noppoint),           &
                 xfoil_options, lift, drag, moment, viscrms)

! Deallocate xfoil variables

  call xfoil_cleanup()

! Deallocate some things

  call deallocate_airfoil(foil)
  deallocate(lift)
  deallocate(drag)
  deallocate(moment)
  deallocate(viscrms)

end program xfoil_only
