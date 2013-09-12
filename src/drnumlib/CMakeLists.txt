SET(drnumlib_SOURCES
    patch.cpp
    blockcfd.cpp
    cartesianpatch.cpp
    codestring.cpp
    timeintegration.cpp
    rungekutta.cpp
    patchgrid.cpp
    patchgroups.cpp
    transformation.cpp
    cartesianraster.cpp
    structuredhexraster.cpp
    raster.cpp
    debug.cpp
    prismaticlayerpatch.cpp
    iteratorfeeder.cpp
    cubeincartisianpatch.cpp
    configmap.cpp
    objectdefinition.cpp
    cartboxobject.cpp
    sphereobject.cpp
    cylinderobject.cpp
    combiobject.cpp
    combiobjector.cpp
    combiobjectand.cpp
    combiobjectandnot.cpp
    blockobject.cpp
    blockobjectbc.cpp
    compressibleeulerbobc.cpp
    coneobject.cpp
    compressiblesolidwallbobc.cpp
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
)

ADD_LIBRARY(drnumlib SHARED ${drnumlib_SOURCES})

SET_TARGET_PROPERTIES(drnumlib PROPERTIES VERSION ${DRNUM_VERSION})
SET_TARGET_PROPERTIES(drnumlib PROPERTIES COMPILE_FLAGS "-Wno-deprecated")

INSTALL(TARGETS drnumlib LIBRARY DESTINATION lib)

