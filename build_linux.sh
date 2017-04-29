#!/bin/bash

INSTALLDIR=$(pwd)/install

if [ -d build ]; then
  rm -rf build
fi
mkdir build

if [ -d install ]; then
  rm -rf install
fi

cd build

  cmake \
    -DCMAKE_INSTALL_PREFIX:PATH="$INSTALLDIR" \
    -DCMAKE_BUILD_TYPE:STRING="Release" \
    ..
  make
  make install

cd ..
