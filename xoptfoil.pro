######################################################################
# QMake project file
######################################################################

TEMPLATE = app
TARGET = xoptfoil
INCLUDEPATH += src/cpp
CPPDIR = src/cpp

# Input
HEADERS += $${CPPDIR}/settingsbrowser.h $${CPPDIR}/optsettings.h $${CPPDIR}/opersettings.h $${CPPDIR}/settingswindow.h $${CPPDIR}/mainwindow.h
SOURCES += $${CPPDIR}/settingsbrowser.cpp $${CPPDIR}/optsettings.cpp $${CPPDIR}/opersettings.cpp $${CPPDIR}/settingswindow.cpp $${CPPDIR}/mainwindow.cpp $${CPPDIR}/gui_main.cpp 
RESOURCES += resources.qrc

QT += widgets
