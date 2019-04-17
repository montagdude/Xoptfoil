# Set up project
cmake_minimum_required(VERSION 2.8.8)
set(XOPTFOIL_VERSION 1.11.1)
project(Xoptfoil)
enable_language(Fortran)

# Some default directories
set(DATADIR ${CMAKE_INSTALL_PREFIX}/share
    CACHE PATH "Directory to install shared data and examples."
)

# Determine compiler
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message("Fortran compiler: " ${CMAKE_Fortran_COMPILER})

# Option to build with/without OpenMP support
set(ENABLE_OPENMP TRUE
    CACHE BOOL "Whether to build with OpenMP support."
)

# Enable OpenMP by default for Release builds, but disable by default for Debug
if (CMAKE_BUILD_TYPE MATCHES "Release")
  set (OPENMP_FLAG "-fopenmp")
  if (NOT ENABLE_OPENMP)
    MESSAGE(WARNING "Disabling OpenMP support since ENABLE_OPENMP=FALSE")
    set (OPENMP_FLAG "")
  endif (NOT ENABLE_OPENMP)
elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
  set (OPENMP_FLAG "")
  if (ENABLE_OPENMP)
    MESSAGE(WARNING "Enabling OpenMP support for debug build since ENABLE_OPENMP=TRUE")
    set (OPENMP_FLAG "-fopenmp")
  endif (ENABLE_OPENMP)
endif (CMAKE_BUILD_TYPE MATCHES "Release")

# Default compiler flags
if (Fortran_COMPILER_NAME MATCHES "gfortran")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${OPENMP_FLAG}")
  set (CMAKE_Fortran_FLAGS_DEBUG "-g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall ${OPENMP_FLAG}")
  set (Fortran_REAL8_FLAG "-fdefault-real-8")
elseif (Fortran_COMPILER_NAME MATCHES "ifort")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${OPENMP_FLAG}")
  set (CMAKE_Fortran_FLAGS_DEBUG "-g -warn all ${OPENMP_FLAG}")
  set (Fortran_REAL8_FLAG "-r8")
else (Fortran_COMPILER_NAME MATCHES "gfortran")
  message(FATAL_ERROR "Fortran compiler not supported.")
endif (Fortran_COMPILER_NAME MATCHES "gfortran")

# Preprocessor definitions
add_definitions(-DPACKAGE_VERSION=\"${XOPTFOIL_VERSION}\")

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
install(FILES src/python/xoptfoil_visualizer.py
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
install(FILES ${AIRFOILS} DESTINATION ${DATADIR}/sample_airfoils)
install(FILES "doc/all_inputs.txt" DESTINATION ${DATADIR})
