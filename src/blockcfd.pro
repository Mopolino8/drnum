TEMPLATE = app
CONFIG += console

INCLUDEPATH += $(VTKINCDIR)
INCLUDEPATH += utility

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
#QMAKE_CXXFLAGS_RELEASE += -g
QMAKE_CXXFLAGS_RELEASE += -O3
QMAKE_CXXFLAGS_RELEASE += -finline-limit=100000
QMAKE_CXXFLAGS_RELEASE += --param large-function-growth=100000
QMAKE_CXXFLAGS_RELEASE += --param inline-unit-growth=100000
QMAKE_CXXFLAGS_RELEASE += -funroll-loops
#QMAKE_CXXFLAGS_RELEASE += -Winline


SOURCES += main.cpp \
    patch.cpp \
    cartesianpatch.cpp \
    reconstruction/limitedreconstruction.cpp \
    blockcfd.cpp \
    timeintegration.cpp \
    patchiterator.cpp \
    rungekutta.cpp \
    patchgrid.cpp \
    utility/List.cc \
    utility/MError.cc \
    math/coordtransform.cpp \
    math/coordtransformvv.cpp \
    transformation.cpp \
    shapes/sphere.cpp \
    shapes/halfspace.cpp

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
    patchgrid.h \
    intercoeff.h \
    intercoeffws.h \
    utility/weightedset.h \
    utility/List.hh \
    utility/MError.hh \
    math/coordtransform.h \
    math/coordtransformvv.h \
    fluxes/ausmdv.h \
    fluxes/ausm.h \
    fluxes/compressiblewallflux.h \
    perfectgas.h \
    iterators/cartesiandirectionalpatchoperation.h \
    shapes/sphere.h \
    fluxes/compressiblefarfieldflux.h \
    reconstruction/immersedboundaryreconstruction.h \
    boundary_conditions/compressibleeulerwall.h \
    shapes/shape.h \
    math/smallsquarematrix.h \
    math/mathvector_structs.h \
    math/mathvector_operators.h \
    math/mathvector_methods.h \
    math/mathvector.h \
    math/linsolve.h \
    transformation.h \
    shapes/halfspace.h \
    compressiblevariables.h
