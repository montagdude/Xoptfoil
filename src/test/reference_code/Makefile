FC=gfortran
OBJ=xfoil_inc.o my_equivalence.o blpar_inc.o xbl_inc.o xfoil_math.o xfoil_geom.o xfoil_bl.o xfoil_aero.o vardef.o math_deps.o xfoil_driver.o parametrization.o memory_util.o airfoil_operations.o

all: test_airfoil

modules: $(OBJ)

test_airfoil: modules test_airfoil.f90
	$(FC) -o test_airfoil test_airfoil.f90 $(OBJ)

clean:
	rm -f *.mod *.o

xfoil_inc.o: xfoil_inc.f90
	$(FC) -c xfoil_inc.f90

my_equivalence.o: my_equivalence.f90
	$(FC) -c my_equivalence.f90

blpar_inc.o: blpar_inc.f90
	$(FC) -c blpar_inc.f90

xbl_inc.o: xbl_inc.f90
	$(FC) -c xbl_inc.f90

xfoil_math.o: xfoil_math.f
	$(FC) -c -fdefault-real-8 xfoil_math.f

xfoil_geom.o: xfoil_geom.f
	$(FC) -c -fdefault-real-8 xfoil_geom.f

xfoil_bl.o: xfoil_bl.f
	$(FC) -c -fdefault-real-8 xfoil_bl.f

xfoil_aero.o: xfoil_aero.f
	$(FC) -c -fdefault-real-8 xfoil_aero.f

vardef.o: vardef.f90
	$(FC) -c vardef.f90

math_deps.o: math_deps.f90
	$(FC) -c math_deps.f90

xfoil_driver.o: xfoil_driver.f90
	$(FC) -c xfoil_driver.f90

parametrization.o: parametrization.f90
	$(FC) -c parametrization.f90

memory_util.o: memory_util.f90
	$(FC) -c memory_util.f90

airfoil_operations.o: airfoil_operations.f90
	$(FC) -c airfoil_operations.f90
