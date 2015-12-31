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

module vardef

! Data structures for airfoil optimization code

  implicit none

  type airfoil_type

    integer :: npoint
    double precision, dimension(:), allocatable :: x, z ! Airfoil coordinates
    double precision :: xle, zle                        ! Leading edge coords
    integer :: leclose                                  ! Index closest to LE
    integer :: addpoint_loc                             ! Whether to add point 
                                                        !  for LE before or
                                                        !  after leclose

  end type airfoil_type

! Global variables (mainly needed to preserve generality of optimization
! routines)

  integer :: noppoint
  integer, parameter :: max_op_points = 30
  double precision, dimension(:), allocatable :: xseedt, xseedb, zseedt, zseedb
  character(7), dimension(max_op_points) :: op_mode
  double precision, dimension(max_op_points) :: op_point, reynolds, mach,      &
                                                flap_degrees, weighting,       &
                                                scale_factor 
  double precision :: x_flap, y_flap
  logical :: use_flap
  character(9), dimension(max_op_points) :: optimization_type

  type(airfoil_type) :: curr_foil
  double precision :: min_thickness, max_thickness, min_moment, min_te_angle,  &
                      growth_allowed
  double precision :: curv_threshold
  integer :: max_curv_reverse
  character(4) :: seed_violation_handling
  character(8) :: moment_constraint_type
  character(11) :: shape_functions
  double precision, dimension(:), allocatable :: xmatcht, xmatchb, zmatcht,    &
                                                 zmatchb
  logical :: match_foils
  logical :: check_curvature
  logical :: symmetrical

  integer :: nparams_top, nparams_bot

!$omp threadprivate(curr_foil)

end module vardef
