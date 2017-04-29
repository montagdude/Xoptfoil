SET INSTALLDIR=%CD%\windows
IF EXIST %INSTALLDIR% RMDIR /S /Q %INSTALLDIR%
IF EXIST build RMDIR /S /Q build
MKDIR %INSTALLDIR%
MKDIR build

CD build
cmake -G "MinGW Makefiles" ^
  -DCMAKE_INSTALL_PREFIX:PATH=%INSTALLDIR% ^
  -DCMAKE_BUILD_TYPE:STRING="Release" ^
  ..

mingw32-make VERBOSE=1
mingw32-make install
CD ..
