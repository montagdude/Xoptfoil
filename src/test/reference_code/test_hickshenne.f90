program test_hickshenne

  use vardef,             only : airfoil_type, xseedt, xseedb, zseedt, zseedb,&
                                 nparams_top, nparams_bot, shape_functions,&
                                 initial_perturb, min_bump_width, curr_foil
  use airfoil_operations, only : get_seed_airfoil, get_split_points,&
                                 split_airfoil
  use memory_util,        only : allocate_airfoil_data, deallocate_airfoil_data
  use parametrization,    only : create_airfoil

  implicit none

  type(airfoil_type) :: seed_airfoil
  double precision :: xoffset, zoffset, foilscale
  character(256) :: airfoil_file
  integer :: i, nargs, pointst, pointsb
  double precision :: t1fact, t2fact
  double precision, dimension(:), allocatable :: modest, modesb
  double precision, dimension(:), allocatable :: zt_new, zb_new

  ! Get airfoil file from command line arguments
  call getarg(1, airfoil_file)

  ! Read seed airfoil
  call get_seed_airfoil("from_file", airfoil_file, seed_airfoil,&
                        xoffset, zoffset, foilscale)

  ! Split seed airfoil into upper and lower surfaces
  call get_split_points(seed_airfoil, pointst, pointsb, .False.)
  allocate(xseedt(pointst))
  allocate(zseedt(pointst))
  allocate(xseedb(pointsb))
  allocate(zseedb(pointsb))
  allocate(zt_new(pointst))
  allocate(zb_new(pointsb))
  call split_airfoil(seed_airfoil, xseedt, xseedb, zseedt, zseedb, .False.)

  ! Set parameters so t1fact and t2fact are 1 in parametrization.f90
  t1fact = 1.d0
  t2fact = 1.d0
  initial_perturb = t1fact*(1.d0 - 0.001d0)
  min_bump_width = 10.d0 - initial_perturb/t2fact

  ! Set up shape functions
  shape_functions = 'hicks-henne'
  nparams_top = 4
  nparams_bot = 4
  call allocate_airfoil_data()
  allocate(modest(nparams_top*3))
  allocate(modesb(nparams_bot*3))
  modest = (/0.05d0, 0.1d0, 0.3d0, -0.03d0, 0.3d0, 0.5d0, 0.1d0, 0.8d0, 1.0d0,&
             0.25d0, 0.5d0, 2.0d0/)
  modesb = (/0.15d0, 0.05d0, 2.0d0, 0.04d0, 0.33d0, 0.3d0, -0.1d0, 0.6d0,&
             0.8d0, -0.03d0, 0.9d0, 3.0d0/)

  ! Create airfoil
  call create_airfoil(xseedt, zseedt, xseedb, zseedb, modest, modesb, zt_new,&
                      zb_new, shape_functions, .false.)

  ! Format coordinates in a single loop in derived type. Also remove translation
  ! and scaling to ensure Cm_x=0.25 doesn't change.
  do i = 1, pointst
    curr_foil%x(i) = xseedt(pointst-i+1)/foilscale - xoffset
    curr_foil%z(i) = zt_new(pointst-i+1)/foilscale - zoffset
  end do
  do i = 1, pointsb-1
    curr_foil%x(i+pointst) = xseedb(i+1)/foilscale - xoffset
    curr_foil%z(i+pointst) = zb_new(i+1)/foilscale - zoffset
  end do

  ! Write coordinates to file
  open(unit=12, file="../data/hickshenne_truthdata.dat", status="replace")
  do i = 1, curr_foil%npoint
    write(12, "(2E21.12)") curr_foil%x(i), curr_foil%z(i)
  end do
  close(12)

  ! Deallocate memory
  deallocate(xseedt, zseedt, xseedb, zseedb, zt_new, zb_new)
  call deallocate_airfoil_data()
  deallocate(modest, modesb)

end program test_hickshenne
