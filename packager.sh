#!/bin/bash

if [ $# -ne 1 ]; then
  echo "Usage: packager.sh xoptfoil_version"
  exit 1
fi

VERSION=$1

SRC_RELEASE="Xoptfoil-${VERSION}-linux_source"
WIN_RELEASE="Xoptfoil-${VERSION}-windows"
rm -rf $SRC_RELEASE
rm -rf $WIN_RELEASE
rm -rf tmp
mkdir $SRC_RELEASE
mkdir $WIN_RELEASE
mkdir tmp

# Top directory: ASCII docs and build scripts
COPYDIR=.
COPYLIST="${COPYDIR}/AUTHORS \
          ${COPYDIR}/CMakeLists.txt \
          ${COPYDIR}/COPYING \
          ${COPYDIR}/ChangeLog \
          ${COPYDIR}/FAQ \
          ${COPYDIR}/INSTALL \
          ${COPYDIR}/README"
cp $COPYLIST tmp

# src directory
COPYDIR='src/fortran src/python'
DESTDIR=src
mkdir tmp/$DESTDIR
cp -r $COPYDIR tmp/$DESTDIR

# doc directory: User guide and input file template
COPYDIR=doc
DESTDIR=doc
COPYLIST="${COPYDIR}/User_Guide.pdf \
          ${COPYDIR}/all_inputs.txt"
mkdir tmp/$DESTDIR
cp $COPYLIST tmp/$DESTDIR

# doc/example_case directory: example cases
COPYDIR=doc/example_case
DESTDIR=doc/example_case
COPYLIST="${COPYDIR}/example_case.pdf \
          ${COPYDIR}/inputs.txt \
          ${COPYDIR}/inputs_withflap.txt"
mkdir tmp/$DESTDIR
cp $COPYLIST tmp/$DESTDIR

# sample_airfoils directory: sample airfoils
COPYDIR=sample_airfoils
cp -r $COPYDIR tmp

# Copy into release directories and remove tmp directory
cp -r tmp/* $SRC_RELEASE
cp -r tmp/* $WIN_RELEASE
rm -rf tmp

# Build scripts
cp build_linux.sh $SRC_RELEASE
cp build_windows.bat $WIN_RELEASE

# Windows executables and DLLs
COPYLIST='windows/bin/* mingw32_dlls/*'
DESTDIR=bin
mkdir $WIN_RELEASE/$DESTDIR
cp $COPYLIST $WIN_RELEASE/$DESTDIR

# Create tarball/zip
rm -f $SRC_RELEASE.tar.gz
rm -f $WIN_RELEASE.zip
tar -czvf $SRC_RELEASE.tar.gz $SRC_RELEASE
zip -rv $WIN_RELEASE.zip $WIN_RELEASE
rm -rf $SRC_RELEASE
rm -rf $WIN_RELEASE
