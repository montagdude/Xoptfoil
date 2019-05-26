program test_airfoil

  use vardef,             only : airfoil_type
  use airfoil_operations, only : get_seed_airfoil

  implicit none

  type(airfoil_type) :: seed_airfoil
  double precision :: xoffset, zoffset, foilscale

  ! Read seed airfoil and write leading edge info
  call get_seed_airfoil("from_file", "../data/mh45_labeled.dat", seed_airfoil,&
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

end program test_airfoil
