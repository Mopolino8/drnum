TEMPLATE = app
CONFIG += console

INCLUDEPATH += $(VTKINCDIR)
INCLUDEPATH += ./utility

LIBS += -L$(VTKLIBDIR)
LIBS += -L/nopt/cuda/5.0/lib64
LIBS += -lQVTK
LIBS += -lvtkCommon
LIBS += -lvtkDICOMParser
LIBS += -lvtkexoIIc
LIBS += -lvtkFiltering
LIBS += -lvtkftgl
LIBS += -lvtkGenericFiltering
LIBS += -lvtkGraphics
LIBS += -lvtkHybrid
LIBS += -lvtkImaging
LIBS += -lvtkIO
LIBS += -lvtkRendering
LIBS += -lvtksys
LIBS += -lvtkVolumeRendering
LIBS += -lvtkWidgets
LIBS += -lgomp


QMAKE_CXXFLAGS += -Wno-deprecated
QMAKE_CXXFLAGS += -fopenmp
#QMAKE_CXXFLAGS += -std=c++0x
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
cuda {
  USEGPU                 = -DGPU
  LIBS                  += -lcudart
  INCLUDEPATH           += .
  INCLUDEPATH           += /usr/include/QtCore
  NVCCFLAGS             += -Xcompiler -fopenmp
  NVCCFLAGS             += -O3
  NVCCFLAGS             += -m64
  NVCCFLAGS             += -arch=sm_20
  NVCCFLAGS             += -Xcompiler -Wno-deprecated
  NVCCFLAGS             += -Xcompiler $$USEGPU -DCUDA
  CUDA_INC               = $$join(INCLUDEPATH,' -I','-I',' ')
  cuda.input             = CUDA_SOURCES
  cuda.output            = ${OBJECTS_DIR}${QMAKE_FILE_BASE}_cuda.o
  cuda.commands          = nvcc -c $$NVCCFLAGS $$CUDA_INC ${QMAKE_FILE_NAME} -o ${QMAKE_FILE_OUT}
  cuda.dependency_type   = TYPE_C
  cuda.depend_command    = nvcc -M $$CUDA_INC $$NVCCFLAGS   ${QMAKE_FILE_NAME}
  QMAKE_EXTRA_COMPILERS += cuda
  QMAKE_CXXFLAGS        += $$USEGPU -DCUDA
  INCLUDEPATH           += $(CUDAINCDIR)
}

cuda_debug {
  USEGPU                 = -DGPU
  LIBS                  += -lcudart
  INCLUDEPATH           += .
  INCLUDEPATH           += /usr/include/QtCore
  NVCCFLAGS             += -Xcompiler -fopenmp
  NVCCFLAGS             += -g -G
  NVCCFLAGS             += -m64
  NVCCFLAGS             += -arch=sm_20
  NVCCFLAGS             += -Xcompiler -Wno-deprecated
  NVCCFLAGS             += -Xcompiler $$USEGPU -DCUDA
  CUDA_INC               = $$join(INCLUDEPATH,' -I','-I',' ')
  cuda.input             = CUDA_SOURCES
  cuda.output            = ${OBJECTS_DIR}${QMAKE_FILE_BASE}_cuda.o
  cuda.commands          = nvcc -c $$NVCCFLAGS $$CUDA_INC ${QMAKE_FILE_NAME} -o ${QMAKE_FILE_OUT}
  cuda.dependency_type   = TYPE_C
  cuda.depend_command    = nvcc -M $$CUDA_INC $$NVCCFLAGS   ${QMAKE_FILE_NAME}
  QMAKE_EXTRA_COMPILERS += cuda
  QMAKE_CXXFLAGS        += $$USEGPU -DCUDA
  INCLUDEPATH           += $(CUDAINCDIR)
}


SOURCES += main.cpp \
    patch.cpp \
    blockcfd.cpp \
    cartesianpatch.cpp \
    codestring.cpp \
    timeintegration.cpp \
    rungekutta.cpp \
    patchgrid.cpp \
    patchgroups.cpp \
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
    rungekuttapg1.cpp \
    prismaticlayerpatch.cpp \
    iteratorfeeder.cpp \
    main.cu

SOURCES -= main.cu

HEADERS += \
    blockcfd.h \
    boundary_conditions/compressibleeulerwall.h \
    cartesianpatch_common.h \
    cartesianpatch.h \
    cartesianraster.h \
    codestring.h \
    compressiblevariables.h \
    debug.h \
    examples/cpuduct.h \
    examples/cpujet.h \
    examples/cpujet_pg.h \
    examples/cpujet_pg_single.h \
    examples/cpujet_pg_dual_1.h \
    examples/cpujet_mb1.h \
    examples/cpujet_mb_grid1.h \
    examples/cpujet_mb_onlygrid.h \
    examples/ffs1.h \
    examples/flatplate.h \
    examples/gpujet.h \
    examples/interpolation1.h \
    examples/jet_common.h \
    examples/jet_pg_common.h \
    examples/jetflux.h \
    examples/kelvinhelmholtz.h \
    examples/overexpandedjet2d.h \
    examples/tjunction.h \
    examples/two_patches_1.h \
    examples/two_patches_2.h \
    examples/two_patches_3.h \
    examples/wedge.h \
    fluxes/ausmdv.h \
    fluxes/ausm.h \
    fluxes/ausmplus.h \
    fluxes/compressiblefarfieldflux.h \
    fluxes/compressibleflux.h \
    fluxes/compressibleslipflux.h \
    fluxes/compressibleviscflux.h \
    fluxes/compressiblewallflux.h \
    fluxes/knp.h \
    fluxes/kt.h \
    fluxes/ktmod.h \
    fluxes/roe.h \
    fluxes/vanleer.h \
    genericoperation.h \
    gpu_cartesianpatch.h \
    gpu_patch.h \
    gpu_rungekutta.h \
    intercoeffpad.h \
    intercoeffws.h \
    iterators/gpu_cartesianiterator.h \
    iterators/gpu_patchiterator.h \
    iterators/patchiterator.h \
    iterators/tpatchiterator.h \
    math/coordtransform.h \
    math/coordtransformvv.h \
    math/linsolve.h \
    math/mathvector.h \
    math/mathvector_methods.h \
    math/mathvector_operators.h \
    math/mathvector_structs.h \
    math/smallsquarematrix.h \
    patch_common.h \
    patchgrid.h \
    patchgroups.h \
    patch.h \
    perfectgas.h \
    raster.h \
    reconstruction/immersedboundaryreconstruction.h \
    reconstruction/minmod.h \
    reconstruction/roelim.h \
    reconstruction/secondorder.h \
    reconstruction/upwind1.h \
    reconstruction/upwind2.h \
    reconstruction/upwindcentral.h \
    reconstruction/vanalbada.h \
    reconstruction/vanleerlim.h \
    rungekutta.h \
    rungekuttapg1.h \
    shapes/box.h \
    shapes/cylindery.h \
    shapes/halfspace.h \
    shapes/noshape.h \
    shapes/shape.h \
    shapes/sphere.h \
    shapes/triangulatedshape.h \
    structuredhexraster.h \
    timeintegration.h \
    tinsecthashraster.h \
    transformation.h \
    utility/List.hh \
    utility/MError.hh \
    utility/namespace_mouse.hh \
    utility/sparseweightedset.h \
    utility/TInsectionList.hh \
    utility/usparseweightedset.h \
    utility/weightedset.h \
    vectorhashraster.h \
    iterators/cartesianiterator.h \
    prismaticlayerpatch.h \
    prismaticlayerpatch_common.h \
    iteratorfeeder.h \
    cudatools.h \
    postprocessingvariables.h \
    stringtools.h \
    donor_t.h \
    external_aero.h


CUDA_SOURCES += main.cu

OTHER_FILES += \
    doc/iterators.dox

