CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# limit configuration types (must be done before project() statement)
SET(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "limited config" FORCE)

PROJECT(DRNUM)

set(A_DEBUG_BUILD OFF CACHE BOOL "Build a debug version?")
if (${A_DEBUG_BUILD})
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
else()
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

if (NOT INITIAL_CONFIG)
  set(INITIAL_CONFIG "..." CACHE INTERNAL "hidden flag...")
  set(CMAKE_INSTALL_PREFIX $ENV{HOME}/local CACHE STRING "Install path prefix, prepended onto install directories." FORCE)
endif()



OPTION(USE_GPU "Use GPU for computation" ON)

SET(CMAKE_VERBOSE_MAKEFILE OFF)

SET(DRNUM_MAJOR_VERSION 0)
SET(DRNUM_MINOR_VERSION 1)
SET(DRNUM_PATCH_VERSION 0)
SET(DRNUM_VERSION ${DRNUM_MAJOR_VERSION}.${DRNUM_MINOR_VERSION}.${DRNUM_PATCH_VERSION})
SET(DRNUM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(DRNUM_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

FIND_PACKAGE(VTK 6.1 NO_MODULE REQUIRED)
#FIND_PACKAGE(VTK 6.0 COMPONENTS vtkRenderingCore vtkChartsCore vtkCommonExecutionModel NO_MODULE REQUIRED)
#include(${VTK_USE_FILE})
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})

set(OPENMP_LIBS, "")
find_package(OpenMP)
if(OpenMP_FOUND)
    add_definitions(-DOPEN_MP)
    set(OPENMP_LIBS, "-lgomp")
endif(OpenMP_FOUND)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  add_definitions(-DMAC_OSX)
endif()

FIND_PACKAGE(MPI REQUIRED)
INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${MPI_C_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})

# Qt5 does not work with (latest) VTK, which links against Qt4
#FIND_PACKAGE(Qt5Core REQUIRED)
#INCLUDE_DIRECTORIES(${Qt5Core_INCLUDE_DIRS})

# Qt4
FIND_PACKAGE(Qt4 COMPONENTS QtCore REQUIRED)
INCLUDE(${QT_USE_FILE})

IF(USE_GPU)
    FIND_PACKAGE(CUDA)
    INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})
    # nvcc flags
    SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} --compiler-options;-fPIC)
    # print kernel statistics
    #LIST(APPEND CUDA_NVCC_FLAGS "-Xptxas -v")
    # generate code for Fermi (2.x) and Kepler (3.x) GPUs
    LIST(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_20,code=compute_20")
    LIST(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_20,code=sm_21")
    LIST(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_30,code=compute_30")
    #LIST(APPEND CUDA_NVCC_FLAGS "-Xptxas -dlcm=ca")
    #LIST(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_35,code=compute_35")
    LIST(APPEND CUDA_NVCC_FLAGS "-Xcompiler -Wno-deprecated")
    LIST(APPEND CUDA_NVCC_FLAGS "-Xcompiler ${OpenMP_CXX_FLAGS}")

    ADD_DEFINITIONS(-DGPU -DCUDA)
ENDIF(USE_GPU)
#ADD_DEFINITIONS(-DDEBUG)

# set include directory, add src directories
INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/drnumlib
    ${CMAKE_CURRENT_BINARY_DIR}/drnumlib
    ${CMAKE_CURRENT_SOURCE_DIR}/shmlib
    ${CMAKE_CURRENT_BINARY_DIR}/shmlib)

ADD_SUBDIRECTORY(shmlib)
ADD_SUBDIRECTORY(drnumlib)
ADD_SUBDIRECTORY(applications)

