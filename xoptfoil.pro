######################################################################
# Project file for qmake
######################################################################

TEMPLATE = app
TARGET = xoptfoil
INCLUDEPATH += src/include

# Input

HEADERS += src/include/optsettings.h \
           src/include/opersettings.h \
           src/include/constrsettings.h \
           src/include/initsettings.h \
           src/include/psosettings.h \
           src/include/gasettings.h \
           src/include/simplexsettings.h \
           src/include/xfanasettings.h \
           src/include/xfpansettings.h \
           src/include/settingsbrowser.h \
           src/include/settingswindow.h \
           src/include/mainwindow.h
SOURCES += src/cpp/optsettings.cpp \
           src/cpp/opersettings.cpp \
           src/cpp/constrsettings.cpp \
           src/cpp/initsettings.cpp \
           src/cpp/psosettings.cpp \
           src/cpp/gasettings.cpp \
           src/cpp/simplexsettings.cpp \
           src/cpp/xfanasettings.cpp \
           src/cpp/xfpansettings.cpp \
           src/cpp/settingsbrowser.cpp \
           src/cpp/settingswindow.cpp \
           src/cpp/mainwindow.cpp \
           src/cpp/gui_main.cpp
RESOURCES += resources.qrc

QT += widgets
