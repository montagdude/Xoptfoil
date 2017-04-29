# Set up project
cmake_minimum_required(VERSION 2.8.8)
set(XOPTFOIL_VERSION 1.11.0)
project(Xoptfoil)
enable_language(Fortran)

# Determine compiler
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message("Fortran compiler: " ${CMAKE_Fortran_COMPILER})

# Default compiler flags
if (Fortran_COMPILER_NAME MATCHES "gfortran")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp")
  set (CMAKE_Fortran_FLAGS_DEBUG "-g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall")
  set (Fortran_REAL8_FLAG "-fdefault-real-8")
elseif (Fortran_COMPILER_NAME MATCHES "ifort")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp")
  set (CMAKE_Fortran_FLAGS_DEBUG "-g -warn all")
  set (Fortran_REAL8_FLAG "-r8")
else (Fortran_COMPILER_NAME MATCHES "gfortran")
  message(FATAL_ERROR "Fortran compiler not supported.")
endif (Fortran_COMPILER_NAME MATCHES "gfortran")

# Some default directories
set(DATADIR ${CMAKE_INSTALL_PREFIX}/share/xoptfoil
    CACHE PATH "Directory to install shared data and examples."
)

# Set source files (exclude programs here)
file(GLOB SOURCESDBLE "src/fortran/*.f90" "src/fortran/*.F90" "src/fortran/xfoil_deps/*.f90")
file(GLOB SOURCESREAL "src/fortran/xfoil_deps/*.f" "src/fortran/naca456_deps/*.f90")
list(REMOVE_ITEM SOURCESDBLE "${CMAKE_SOURCE_DIR}/src/fortran/main.f90")
list(REMOVE_ITEM SOURCESDBLE "${CMAKE_SOURCE_DIR}/src/fortran/xfoil_only.f90")

# Real-8 flag for xfoil and naca456 source files
set_source_files_properties(${SOURCESREAL} PROPERTIES COMPILE_FLAGS ${Fortran_REAL8_FLAG})

# CMake "object libraries" which allow sharing object files for multiple executables
# https://cmake.org/Wiki/CMake/Tutorials/Object_Library
# (Introduced in CMake 2.8.8)
add_library(MYOBJLIB OBJECT ${SOURCESDBLE} ${SOURCESREAL})

# Build executables
add_executable("xoptfoil" $<TARGET_OBJECTS:MYOBJLIB> "src/fortran/main.f90")
add_executable("xfoil_only" $<TARGET_OBJECTS:MYOBJLIB> "src/fortran/xfoil_only.f90")

# Installation
install(TARGETS xoptfoil xfoil_only DESTINATION bin)
install(FILES xoptfoil_visualizer
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ             GROUP_EXECUTE
                    WORLD_READ             WORLD_EXECUTE
        DESTINATION bin)
install(FILES "doc/User_Guide.pdf" DESTINATION ${DATADIR})
install(FILES "doc/example_case/example_case.pdf"
              "doc/example_case/inputs.txt" 
              "doc/example_case/inputs_withflap.txt" 
              DESTINATION ${DATADIR}/example_case)
file(GLOB AIRFOILS "sample_airfoils/*.dat")
install(FILES ${AIRFOILS} DESTINATION ${DATADIR})
install(FILES "inputs.txt" DESTINATION ${DATADIR})
