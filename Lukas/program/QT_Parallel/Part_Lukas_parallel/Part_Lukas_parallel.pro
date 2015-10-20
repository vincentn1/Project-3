TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -fopenmp
QMAKE_CXXFLAGS += -O3

LIBS += -fopenmp

SOURCES += main.cpp

HEADERS +=

