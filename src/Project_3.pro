TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -fopenmp -larmadillo -llapack -lblas

QMAKE_CXXFLAGS_RELEASE += -DARMA_NO_DEBUG

SOURCES += main.cpp \
    lib.cpp \
    System.cpp \
    Potential.cpp \
    BasisFunctions.cpp \
    AlphaHarmonicOscillator.cpp \
    Jastrow.cpp \
    GaussianDeviate.cpp

HEADERS += \
    lib.h \
    System.h \
    Potential.h \
    Structs.h \
    BasisFunctions.h \
    AlphaHarmonicOscillator.h \
    Jastrow.h \
    GaussianDeviate.h

OTHER_FILES += \
    lib.o \
    _lib.so \
    lib.i

QMAKE_CXXFLAGS += -std=c++11 -fopenmp
