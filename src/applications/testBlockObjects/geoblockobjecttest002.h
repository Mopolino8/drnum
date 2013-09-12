// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of DrNUM.                                          +
// +                                                                      +
// + Copyright 2013 numrax GmbH, enGits GmbH                              +
// +                                                                      +
// + DrNUM is free software: you can redistribute it and/or modify        +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + DrNUM is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "blockcfd.h"

// Includes to build a geometric test case

// some spheres
SphereObject obj_sph_1;
SphereObject obj_sph_2;
SphereObject obj_sph_3;
obj_sph_1.setParams(4., 4., 4., 2.5);
obj_sph_2.setParams(6., 6., 6., 2.5);
obj_sph_3.setParams(4., 5., 6., 2.5);
// some cylinders
CylinderObject obj_cyl_1;
obj_cyl_1.setParams(1., 5., 5.,
                    6., 0., 0.,
                    1.);

CombiObjectAnd object_inter1(&obj_sph_1, &obj_sph_2);
CombiObjectAnd object_inter2(&obj_sph_2, &obj_sph_3);
CombiObjectAnd object_inter3(&obj_sph_3, &obj_sph_1);
CombiObjectOr object_inter4(&object_inter1, &object_inter2);
CombiObjectOr object_inter5(&object_inter4, &object_inter3);
CombiObjectAndNot object(&object_inter5, &obj_cyl_1);

gridfile = "patches/standard.grid";
