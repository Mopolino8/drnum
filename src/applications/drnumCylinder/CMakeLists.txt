SET(drnumCylinder_CC_SOURCES main.cpp)
SET(DRNUM_USED_LIBS drnumlib shmlib)

IF(USE_GPU)
    SET(drnumCylinder_CU_SOURCES main.cu)
    cuda_add_library(drnumCylinder_cuda ${drnumCylinder_CU_SOURCES} SHARED) # OPTIONS -G)
    SET(DRNUM_USED_LIBS drnumCylinder_cuda ${DRNUM_USED_LIBS})
ENDIF(USE_GPU)


ADD_EXECUTABLE(drnumCylinder ${drnumCylinder_CC_SOURCES})
ADD_DEPENDENCIES(drnumCylinder ${DRNUM_USED_LIBS})
TARGET_LINK_LIBRARIES(drnumCylinder ${DRNUM_USED_LIBS} ${QT_LIBRARIES} ${VTK_LIBRARIES} ${MPI_LIBRARIES} "-lgomp")
#QT5_USE_MODULES(drnumBasicAero Core)

SET_TARGET_PROPERTIES(drnumCylinder
    PROPERTIES
    LINKER_LANGUAGE CXX
    PREFIX "")

SET_TARGET_PROPERTIES(drnumCylinder
    PROPERTIES
    VERSION ${DRNUM_VERSION})

INSTALL(TARGETS drnumCylinder RUNTIME DESTINATION bin)
IF(USE_GPU)
    INSTALL(TARGETS drnumCylinder_cuda LIBRARY DESTINATION lib)
ENDIF(USE_GPU)

