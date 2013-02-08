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
QMAKE_CXXFLAGS += -fopenmp
#QMAKE_CXXFLAGS_RELEASE += -g
QMAKE_CXXFLAGS_RELEASE += -O3
QMAKE_CXXFLAGS_RELEASE += -finline-limit=100000
QMAKE_CXXFLAGS_RELEASE += --param large-function-growth=100000
QMAKE_CXXFLAGS_RELEASE += --param inline-unit-growth=100000
#QMAKE_CXXFLAGS_RELEASE += -funroll-loops
#QMAKE_CXXFLAGS_RELEASE += -Winline


# CUDA
# ====
#
LIBS                  += -lcudart
LIBS                  += -lgomp
INCLUDEPATH           += /opt/CUDA-4.2
NVCCFLAGS             += -Xcompiler -fopenmp -O3 -m64 -arch=sm_13# -Xptxas -dlcm=cg #NVCCFLAGS  += -g -G
CUDA_INC               = $$join(INCLUDEPATH,' -I','-I',' ')
cuda.input             = CUDA_SOURCES
cuda.output            = ${OBJECTS_DIR}${QMAKE_FILE_BASE}_cuda.o
cuda.commands          = nvcc -c $$NVCCFLAGS $$CUDA_INC ${QMAKE_FILE_NAME} -o ${QMAKE_FILE_OUT}
cuda.dependency_type   = TYPE_C
cuda.depend_command    = nvcc -M $$CUDA_INC $$NVCCFLAGS   ${QMAKE_FILE_NAME}
QMAKE_EXTRA_COMPILERS += cuda



SOURCES += main.cpp \
    patch.cpp \
    cartesianpatch.cpp \
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
    shapes/halfspace.cpp \
    cartesianraster.cpp \
    utility/TInsectionList.cc \
    structuredhexraster.cpp \
    raster.cpp \
    shapes/triangulatedshape.cpp \
    shapes/box.cpp \
    debug.cpp \
    shapes/cylindery.cpp \
    rungekuttapg1.cpp

HEADERS += \
    patch.h \
    cartesianpatch.h \
    blockcfd.h \
    reconstruction/upwind1.h \
    fluxes/ausmplus.h \
    fluxes/ausmplus.h \
    reconstruction/upwind2.h \
    iterators/cartesianpatchoperation.h \
    timeintegration.h \
    patchiterator.h \
    iterators/cartesianstandarditerator.h \
    iterators/cartesianstandardpatchoperation.h \
    rungekutta.h \
    patchgrid.h \
    intercoeffpad.h \
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
    compressiblevariables.h \
    cartesianraster.h \
    utility/TInsectionList.hh \
    vectorhashraster.h \
    tinsecthashraster.h \
    utility/sparseweightedset.h \
    utility/usparseweightedset.h \
    structuredhexraster.h \
    raster.h \
    shapes/triangulatedshape.h \
    shapes/noshape.h \
    shapes/box.h \
    reconstruction/minmod.h \
    reconstruction/vanalbada.h \
    fluxes/compressibleflux.h \
    fluxes/kt.h \
    fluxes/knp.h \
    reconstruction/upwindcentral.h \
    reconstruction/vanleerlim.h \
    debug.h \
    fluxes/roe.h \
    fluxes/vanleer.h \
    reconstruction/roelim.h \
    fluxes/compressibleviscflux.h \
    shapes/cylindery.h \
    reconstruction/secondorder.h \
    examples/ffs1.h \
    examples/interpolation1.h \
    examples/overexpandedjet2d.h \
    examples/kelvinhelmholtz.h \
    examples/flatplate.h \
    genericoperation.h \
    fluxes/ktmod.h \
    examples/wedge.h \
    examples/tjunction.h \
    examples/two_patches_1.h \
    examples/two_patches_2.h \
    examples/two_patches_3.h \
    rungekuttapg1.h \
    gpu_cartesianpatch.h \
    cartesianpatch_common.h \
    patch_common.h \
    blockcfd_cuda.h \
    examples/cpuduct.h \
    examples/cpujet.h \
    gpu_patch.h
