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
mkdir $dir/src
mkdir $dir/bin
mkdir $dir/doc
mkdir $dir/doc/example_case
cp -r $cdir/src/fortran/ $dir/src
cp -r $cdir/license/ $dir
cp -r $cdir/sample_airfoils/ $dir
cp $cdir/bin/inputs.txt $dir/bin
cp $cdir/bin/Makefile_* $dir/bin
cp $cdir/bin/windows_cmd_prompt_here.bat $dir/bin
cp $cdir/bin/design_visualizer.py $dir/bin
cp -r $cdir/bin/x86-64/ $dir/bin
cp $cdir/README $dir
cp $cdir/INSTALL $dir
cp $cdir/CHANGELOG $dir
cp $cdir/FAQ $dir
cp $cdir/doc/User_Guide.pdf $dir/doc
cp $cdir/doc/example_case/example_case.pdf $dir/doc/example_case
cp $cdir/doc/example_case/inputs*.txt $dir/doc/example_case

# Zip it
zip -r ${dir}.zip $dir

# Remove $dir once it's been zipped
rm -rf $dir
