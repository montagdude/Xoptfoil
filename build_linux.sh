#!/bin/bash

INSTALLDIR=$(pwd)/linux

rm -rf build
mkdir build
rm -rf linux

cd build

  cmake \
    -DCMAKE_INSTALL_PREFIX:PATH="$INSTALLDIR" \
    -DCMAKE_BUILD_TYPE:STRING="Release" \
    ..
  make VERBOSE=1 || exit 1
  make install || exit 1

cd ..

