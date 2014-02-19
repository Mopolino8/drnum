SET(drnumBaseFlowDemo_CC_SOURCES main.cpp baseflowflux.h)
SET(DRNUM_USED_LIBS drnumlib shmlib)

IF(USE_GPU)
    SET(drnumBaseFlowDemo_CU_SOURCES main.cu)
    cuda_add_library(drnumBaseFlowDemo_cuda ${drnumBaseFlowDemo_CU_SOURCES} SHARED) # OPTIONS -G)
    SET(DRNUM_USED_LIBS drnumBaseFlowDemo_cuda ${DRNUM_USED_LIBS})
ENDIF(USE_GPU)


ADD_EXECUTABLE(drnumBaseFlowDemo ${drnumBaseFlowDemo_CC_SOURCES})
ADD_DEPENDENCIES(drnumBaseFlowDemo ${DRNUM_USED_LIBS})
TARGET_LINK_LIBRARIES(drnumBaseFlowDemo ${DRNUM_USED_LIBS} ${QT_LIBRARIES} ${VTK_LIBRARIES} ${MPI_LIBRARIES} "-lgomp")
#QT5_USE_MODULES(drnumJetDemo Core)

SET_TARGET_PROPERTIES(drnumBaseFlowDemo
    PROPERTIES
    LINKER_LANGUAGE CXX
    PREFIX "")

SET_TARGET_PROPERTIES(drnumBaseFlowDemo
    PROPERTIES
    VERSION ${DRNUM_VERSION})

INSTALL(TARGETS drnumBaseFlowDemo RUNTIME DESTINATION bin)
IF(USE_GPU)
    INSTALL(TARGETS drnumBaseFlowDemo_cuda LIBRARY DESTINATION lib)
ENDIF(USE_GPU)

