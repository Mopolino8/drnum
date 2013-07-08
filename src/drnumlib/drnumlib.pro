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
TEMPLATE = lib
drnumlib.path = ../../lib
drnumlib.files = libdrnumlib.*
INSTALLS += drnumlib

include (../drnum.pri)

SOURCES += \
    patch.cpp \
    blockcfd.cpp \
    cartesianpatch.cpp \
    codestring.cpp \
    timeintegration.cpp \
    rungekutta.cpp \
    patchgrid.cpp \
    patchgroups.cpp \
    math/coordtransform.cpp \
    math/coordtransformvv.cpp \
    transformation.cpp \
    cartesianraster.cpp \
    structuredhexraster.cpp \
    raster.cpp \
    debug.cpp \
    prismaticlayerpatch.cpp \
    iteratorfeeder.cpp \
    blockobject.cpp \
    sphereobject.cpp \
    cartboxobject.cpp \
    cubeincartisianpatch.cpp \
    configmap.cpp


HEADERS += \
    blockcfd.h \
    cartesianpatch_common.h \
    cartesianpatch.h \
    cartesianraster.h \
    codestring.h \
    compressiblevariables.h \
    debug.h \
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
    structuredhexraster.h \
    timeintegration.h \
    tinsecthashraster.h \
    transformation.h \
    sparseweightedset.h \
    usparseweightedset.h \
    weightedset.h \
    vectorhashraster.h \
    iterators/cartesianiterator.h \
    prismaticlayerpatch.h \
    iteratorfeeder.h \
    cudatools.h \
    postprocessingvariables.h \
    stringtools.h \
    donor_t.h \
    cubeincartisianpatch.h \
    configmap.h \
    blockobject.h \
    sphereobject.h \
    cartboxobject.h \
    cubeincartisianpatch.h

