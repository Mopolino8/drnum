# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of DrNUM.                                          +
# +                                                                      +
# + Copyright 2013 numrax GmbH, enGits GmbH                              +
# +                                                                      +
# + DrNUM is free software: you can redistribute it and/or modify        +
# + it under the terms of the GNU General Public License as published by +
# + the Free Software Foundation, either version 3 of the License, or    +
# + (at your option) any later version.                                  +
# +                                                                      +
# + enGrid is distributed in the hope that it will be useful,            +
# + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
# + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
# + GNU General Public License for more details.                         +
# +                                                                      +
# + You should have received a copy of the GNU General Public License    +
# + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
# +                                                                      +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
include (../../drnum.pri)

LIBS += -L../../drnumlib
LIBS += -ldrnumlib

LIBS += -L$(VTKLIBDIR)
LIBS += -L$(CUDALIBDIR)
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

INCLUDEPATH += ../../drnumlib

# CUDA
# ====
#
cuda {
  USEGPU                 = -DGPU
  LIBS                  += -L$(CUDALIBDIR)
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
  LIBS                  += -L$(CUDALIBDIR)
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

