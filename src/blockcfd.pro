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

SOURCES += main.cpp \
    patch.cpp \
    cartesianpatch.cpp \
    fluxes/compressibleflux.cpp \
    reconstruction/limitedreconstruction.cpp

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
    fluxes/ausm.h \
    reconstruction/upwind2.h \
    reconstruction/limitedreconstruction.h

