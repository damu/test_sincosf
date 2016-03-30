TEMPLATE = app
CONFIG += console C++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp
QMAKE_CXXFLAGS_DEBUG += -g
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS += -O3 -ffast-math

include(deployment.pri)
qtcAddDeployment()

