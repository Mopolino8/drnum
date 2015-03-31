include(${VTK_USE_FILE})

SET(drnumlib_SOURCES
    blockcfd.cpp
    blockobjectbc.cpp
    blockobject.cpp
    cartboxobject.cpp
    cartesiancycliccopy.h
    cartesianpatch.cpp
    cartesianpatch_common.h
    cartesianraster.cpp
    codestring.cpp
    combiobjectand.cpp
    combiobjectandnot.cpp
    combiobject.cpp
    combiobjector.cpp
    compressibleeulerbobc.cpp
    compressiblesolidwallbobc.cpp
    coneobject.cpp
    configmap.cpp
    containertricks.h
    cubeincartisianpatch.cpp
    cylinderincartisianpatch.cpp
    cudatools.h
    cylinderobject.cpp
    debug.cpp
    discretelevelset.h
    donor_t.h
    drnum.h
    externalexchangelist.cpp
    geometrytools.cpp
    gpu_cartesianpatch.h
    gpu_cartesianlevelsetbc.h
    gpu_cylinderincartisianpatch.h
    gpu_patch.h
    iteratorfeeder.cpp
    levelsetdefinition.cpp
    levelsetobject.cpp
    levelsetobjectbc.cpp
    mpicommunicator.cpp
    objectdefinition.cpp
    patch.cpp
    patch_common.h
    patchgrid.cpp
    patchgroups.cpp
    perfectgas.h
    prismaticlayerpatch.cpp
    raster.cpp
    rungekutta.cpp
    sphereobject.cpp
    spherelevelset.cpp
    splitface_t.h
    stringtools.h
    structuredhexraster.cpp
    utilities.cpp
    timeintegration.cpp
    transformation.cpp
    iterators/cartesianiterator.h
    iterators/gpu_cartesianiterator.h
    iterators/gpu_patchiterator.h
    iterators/patchiterator.h
    iterators/tpatchiterator.h
    fluxes/ausmdv.h
    fluxes/ausm.h
    fluxes/ausmplus.h
    fluxes/compressiblefarfieldflux.h
    fluxes/compressibleflux.h
    fluxes/compressibleslipflux.h
    fluxes/compressiblesmagorinskyflux.h
    fluxes/compressibleviscflux.h
    fluxes/compressiblewallflux.h
    fluxes/knp.h
    fluxes/kt.h
    fluxes/ktmod.h
    fluxes/roe.h
    fluxes/vanleer.h
    math/coordtransform.cpp
    math/coordtransformvv.cpp
    math/linsolve.h
    math/mathvector.h
    math/mathvector_methods.h
    math/mathvector_operators.h
    math/mathvector_structs.h
    math/smallsquarematrix.h
    reconstruction/minmod.h
    reconstruction/roelim.h
    reconstruction/secondorder.h
    reconstruction/upwind1.h
    reconstruction/upwind2.h
    reconstruction/upwindcentral.h
    reconstruction/vanalbada.h
    reconstruction/vanleerlim.h
    cpu_levelsetobjectbc.h
    lsobccompressibleswallop.h
    timeaverage.h
    timeaverage.cpp
)

ADD_LIBRARY(drnumlib SHARED ${drnumlib_SOURCES})

SET_TARGET_PROPERTIES(drnumlib PROPERTIES VERSION ${DRNUM_VERSION})
SET_TARGET_PROPERTIES(drnumlib PROPERTIES COMPILE_FLAGS "-Wno-deprecated")

TARGET_LINK_LIBRARIES(drnumlib ${DRNUM_USED_LIBS} ${QT_LIBRARIES} ${VTK_LIBRARIES} ${MPI_LIBRARIES} ${OPENMP_LIBS} shmlib)

INSTALL(TARGETS drnumlib LIBRARY DESTINATION lib)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET_TARGET_PROPERTIES(drnumlib PROPERTIES MACOSX_RPATH 1)
endif()


