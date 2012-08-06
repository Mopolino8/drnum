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



SOURCES += main.cpp \
    patch.cpp \
    cartesianpatch.cpp \
    reconstruction/limitedreconstruction.cpp \
    blockcfd.cpp \
    compressiblecartesianpatch.cpp \
    timeintegration.cpp \
    patchiterator.cpp \
    rungekutta.cpp \
    patchgrid.cpp \
    utility/List.cc \
    utility/MError.cc \
    math/coordtransform.cpp \
    math/coordtransformvv.cpp

HEADERS += \
    patch.h \
    cartesianpatch.h \
    blockcfd.h \
    compressiblecartesianpatch.h \
    ausmtools.h \
    reconstruction/upwind1.h \
    fluxes/compressibleflux.h \
    fluxes/ausmbase.h \
    fluxes/ausmplus.h \
    fluxes/ausmplus.h \
    compressibleobject.h \
    reconstruction/upwind2.h \
    reconstruction/limitedreconstruction.h \
    cartesianpatchoperation.h \
    compressiblecartesianpatchoperation.h \
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
    math/coordtransformvv.h
