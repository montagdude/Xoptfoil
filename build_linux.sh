#!/bin/bash

INSTALLDIR=$(pwd)/linux

if [ -d build ]; then
  rm -rf build
fi
mkdir build

if [ -d linux ]; then
  rm -rf linux
fi

cd build

  cmake \
    -DCMAKE_INSTALL_PREFIX:PATH="$INSTALLDIR" \
    -DCMAKE_BUILD_TYPE:STRING="Release" \
    ..
  make VERBOSE=1 || exit 1
  make install || exit 1

cd ..
