#!/bin/bash

# Check input
version=$1
if [ "$version" == "" ]; then
  echo
  echo "Usage: ./packager.sh {xoptfoil_version}"
  echo
  exit
fi

# Create directory
dir="Xoptfoil_${version}"
cdir="."
if [ -d "$dir" ]; then
  rm -rf $dir
fi

# Copy files
mkdir $dir
mkdir -p $dir/src/fortran
mkdir $dir/bin
mkdir $dir/doc
mkdir $dir/doc/example_case
mkdir -p $dir/references/genetic_algorithm $dir/references/parametrization $dir/references/simplex_search
cp $cdir/src/fortran/*.f90 $dir/src/fortran
cp $cdir/src/fortran/*.F90 $dir/src/fortran
cp -r $cdir/src/fortran/xfoil_deps/ $dir/src/fortran
cp -r $cdir/license/ $dir
cp -r $cdir/sample_airfoils/ $dir
cp $cdir/bin/inputs.txt $dir/bin
cp $cdir/bin/inputs_xfoil_only.txt $dir/bin
cp $cdir/bin/Makefile_* $dir/bin
cp $cdir/bin/windows_cmd_prompt_here.bat $dir/bin
cp $cdir/bin/design_visualizer.py $dir/bin
cp -r $cdir/bin/x86-64/ $dir/bin
cp $cdir/README $dir
cp $cdir/INSTALL $dir
cp $cdir/doc/*.pdf $dir/doc
cp $cdir/doc/example_case/example_case.pdf $dir/doc/example_case
cp $cdir/doc/example_case/inputs*.txt $dir/doc/example_case
cp $cdir/references/genetic_algorithm/GAs.odp $dir/references/genetic_algorithm
cp $cdir/references/parametrization/* $dir/references/parametrization
cp $cdir/references/simplex_search/* $dir/references/simplex_search

# Zip it
zip -r ${dir}.zip $dir

# Remove $dir once it's been zipped
rm -rf $dir
