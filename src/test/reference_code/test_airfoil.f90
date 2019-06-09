program test_airfoil

  use vardef,             only : airfoil_type, xseedt, xseedb, zseedt, zseedb
  use airfoil_operations, only : get_seed_airfoil, get_split_points,&
                                 split_airfoil

  implicit none

  type(airfoil_type) :: seed_airfoil
  double precision :: xoffset, zoffset, foilscale
  character(256) :: airfoil_file
  integer :: i, nargs, pointst, pointsb

  ! Get airfoil file from command line arguments
  call getarg(1, airfoil_file)

  ! Read seed airfoil and write leading edge info
  call get_seed_airfoil("from_file", airfoil_file, seed_airfoil,&
                        xoffset, zoffset, foilscale)
  open(unit=12, file="../data/le_truthdata.dat", status="replace")
  write(12, "(E15.8)") seed_airfoil%xle
  write(12, "(E15.8)") seed_airfoil%zle
  write(12, "(I8)") seed_airfoil%leclose
  write(12, "(I8)") seed_airfoil%addpoint_loc
  write(12, "(ES15.8)") xoffset
  write(12, "(ES15.8)") zoffset
  write(12, "(ES15.8)") foilscale
  close(12)

  ! Split seed airfoil into upper and lower surfaces
  call get_split_points(seed_airfoil, pointst, pointsb, .False.)
  allocate(xseedt(pointst))
  allocate(zseedt(pointst))
  allocate(xseedb(pointsb))
  allocate(zseedb(pointsb))
  call split_airfoil(seed_airfoil, xseedt, xseedb, zseedt, zseedb, .False.)

  ! Write upper and lower surfaces to file
  open(unit=13, file="../data/splitfoil_truthdata.dat", status="replace")
  write(13, "(I8)") pointst
  do i = 1, pointst
    write(13, "(2E21.12)") xseedt(i), zseedt(i)
  end do
  write(13, "(I8)") pointsb
  do i = 1, pointsb
    write(13, "(2E21.12)") xseedb(i), zseedb(i)
  end do
  close(13)

  deallocate(xseedt, zseedt, xseedb, zseedb)

end program test_airfoil
