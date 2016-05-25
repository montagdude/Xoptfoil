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
  character(8), dimension(max_op_points) :: flap_selection
  double precision, dimension(max_op_points) :: op_point, reynolds, mach,      &
                                           flap_degrees, weighting, scale_factor 
  double precision :: x_flap, y_flap
  logical :: use_flap
  character(9), dimension(max_op_points) :: optimization_type
  integer :: nflap_optimize          ! Number of operating points where flap 
                                     !   setting will be optimized
  integer, dimension(max_op_points) :: flap_optimize_points

  type(airfoil_type) :: curr_foil
  double precision :: min_thickness, max_thickness, min_te_angle,              &
                      growth_allowed, min_flap_degrees, max_flap_degrees,      &
                      min_camber, max_camber
  double precision :: curv_threshold
  integer :: max_curv_reverse_top, max_curv_reverse_bot
  character(4) :: seed_violation_handling
  character(8), dimension(max_op_points) :: moment_constraint_type
  double precision, dimension(max_op_points) :: min_moment
  character(11) :: shape_functions
  double precision, dimension(:), allocatable :: xmatcht, xmatchb, zmatcht,    &
                                                 zmatchb
  logical :: match_foils
  logical :: check_curvature
  logical :: symmetrical

  integer :: nparams_top, nparams_bot
  double precision :: initial_perturb
  double precision :: min_bump_width

  character(80) :: output_prefix

!$omp threadprivate(curr_foil)

end module vardef
