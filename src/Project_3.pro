TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    lib.cpp \
    System.cpp \
    Potential.cpp

HEADERS += \
    lib.h \
    System.h \
    Potential.h \
    Structs.h

OTHER_FILES += \
    lib.o \
    _lib.so \
    lib.i

QMAKE_CXXFLAGS += -std=c++11
