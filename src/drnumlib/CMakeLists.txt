SET(drnumlib_SOURCES
    patch.cpp
    blockcfd.cpp
    cartesianpatch.cpp
    codestring.cpp
    timeintegration.cpp
    rungekutta.cpp
    patchgrid.cpp
    patchgroups.cpp
    math/coordtransform.cpp
    math/coordtransformvv.cpp
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
    )

ADD_LIBRARY(drnumlib SHARED ${drnumlib_SOURCES})

SET_TARGET_PROPERTIES(drnumlib
    PROPERTIES
    VERSION ${DRNUM_VERSION})

INSTALL(TARGETS drnumlib LIBRARY DESTINATION lib)

