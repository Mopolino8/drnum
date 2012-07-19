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
    compressiblecartesianpatch.cpp

HEADERS += \
    patch.h \
    cartesianpatch.h \
    blockcfd.h \
    compressiblecartesianpatch.h \
    fluxes/ausm_plus_x.h \
    ausmtools.h \
    reconstruction/minmod.h \
    reconstruction/first_order.h \
    reconstruction/second_order.h \
    reconstruction/van_albada1.h \
    reconstruction/van_albada2.h \
    reconstruction/van_leer.h \
    reconstruction/superbee.h

