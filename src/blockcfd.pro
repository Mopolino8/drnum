TEMPLATE = app
CONFIG += console

INCLUDEPATH += $(VTKINCDIR)

LIBS        += -L$(VTKLIBDIR)
LIBS        += -lQVTK
LIBS        += -lvtkCommon
LIBS        += -lvtkDICOMParser
LIBS        += -lvtkexoIIc
LIBS        += -lvtkFiltering
LIBS        += -lvtkftgl
LIBS        += -lvtkGenericFiltering
LIBS        += -lvtkGraphics
LIBS        += -lvtkHybrid
LIBS        += -lvtkImaging
LIBS        += -lvtkIO
LIBS        += -lvtkRendering
LIBS        += -lvtksys
LIBS        += -lvtkVolumeRendering
LIBS        += -lvtkWidgets

QMAKE_CXXFLAGS += -Wno-deprecated
QMAKE_CXXFLAGS_RELEASE += -g
QMAKE_CXXFLAGS_RELEASE += -O3
QMAKE_CXXFLAGS_RELEASE += -finline-limit=100000
QMAKE_CXXFLAGS_RELEASE += --param large-function-growth=100000
QMAKE_CXXFLAGS_RELEASE += --param inline-unit-growth=100000
QMAKE_CXXFLAGS_RELEASE += -funroll-loops
#QMAKE_CXXFLAGS_RELEASE += -Winline

QMAKE_LFLAGS_RELEASE += -O3

SOURCES += main.cpp \
    patch.cpp \
    cartesianpatch.cpp \
    reconstruction/limitedreconstruction.cpp \
    blockcfd.cpp \
    timeintegration.cpp \
    patchiterator.cpp \
    rungekutta.cpp

HEADERS += \
    patch.h \
    cartesianpatch.h \
    blockcfd.h \
    reconstruction/upwind1.h \
    fluxes/ausmbase.h \
    fluxes/ausmplus.h \
    fluxes/ausmplus.h \
    reconstruction/upwind2.h \
    reconstruction/limitedreconstruction.h \
    iterators/cartesianpatchoperation.h \
    timeintegration.h \
    patchiterator.h \
    iterators/cartesianstandarditerator.h \
    iterators/cartesianstandardpatchoperation.h \
    rungekutta.h \
    fluxes/ausmdv.h \
    fluxes/ausm.h \
    fluxes/compressiblewallflux.h \
    perfectgas.h \
    iterators/cartesiandirectionalpatchoperation.h \
    shapes/sphere.h \
    fluxes/compressiblefarfieldflux.h \
    reconstruction/immersedboundaryreconstruction.h \
    boundary_conditions/compressibleeulerwall.h \
    shapes/shape.h

