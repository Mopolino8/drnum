CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# limit configuration types (must be done before project() statement)
SET(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "limited config" FORCE)

PROJECT(DRNUM)
OPTION(USE_GPU "Use GPU for computation" ON)

SET(CMAKE_VERBOSE_MAKEFILE off)

SET(DRNUM_MAJOR_VERSION 0)
SET(DRNUM_MINOR_VERSION 1)
SET(DRNUM_PATCH_VERSION 0)
SET(DRNUM_VERSION ${DRNUM_MAJOR_VERSION}.${DRNUM_MINOR_VERSION}.${DRNUM_PATCH_VERSION})
SET(DRNUM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(DRNUM_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

FIND_PACKAGE(VTK COMPONENTS vtkChartsCore NO_MODULE REQUIRED)
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})

# Qt5 does not work with (latest) VTK, which links against Qt4
#FIND_PACKAGE(Qt5Core REQUIRED)
#INCLUDE_DIRECTORIES(${Qt5Core_INCLUDE_DIRS})

# Qt4
FIND_PACKAGE(Qt4 COMPONENTS QtCore REQUIRED)
INCLUDE(${QT_USE_FILE})

IF(USE_GPU)
    FIND_PACKAGE(CUDA REQUIRED)
    INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})

    # nvcc flags
    SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} --compiler-options;-fPIC)
    # print kernel statistics
    LIST(APPEND CUDA_NVCC_FLAGS "-Xptxas -v")
    # generate code for Fermi (2.x) and Kepler (3.x) GPUs
    LIST(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_20,code=compute_20")
    LIST(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_20,code=sm_21")
    LIST(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_30,code=compute_30")
    #LIST(APPEND CUDA_NVCC_FLAGS "-gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_35,code=compute_35")

    ADD_DEFINITIONS(-DGPU -DCUDA)
ENDIF(USE_GPU)
ADD_DEFINITIONS(-DDEBUG)

# set include directory, add src directories
INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/drnumlib
    ${CMAKE_CURRENT_BINARY_DIR}/drnumlib)

ADD_SUBDIRECTORY(applications)
ADD_SUBDIRECTORY(drnumlib)
ADD_SUBDIRECTORY(shmlib)
