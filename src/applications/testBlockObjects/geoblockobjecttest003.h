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

// Model a desk with a vase on it
CartboxObject plate;
CartboxObject leg_1;
CartboxObject leg_2;
CartboxObject leg_3;
CartboxObject leg_4;
SphereObject vase_body;
CylinderObject vase_inner;

real x_p_min = 2.0;
real x_p_max = 7.0;
real y_p_min = 3.0;
real y_p_max = 7.0;
real leg_thick = 0.5;

plate.setParams(x_p_min, x_p_max,
                y_p_min, y_p_max,
                2.5, 3.0);

leg_1.setParams(x_p_min,             x_p_min + leg_thick,
                y_p_min,             y_p_min + leg_thick,
                -1.0, 2.75);

leg_2.setParams(x_p_min,             x_p_min + leg_thick,
                y_p_max - leg_thick, y_p_max,
                -1.0, 2.75);

leg_3.setParams(x_p_max - leg_thick, x_p_max,
                y_p_min,             y_p_min + leg_thick,
                -1.0, 2.75);

leg_4.setParams(x_p_max - leg_thick, x_p_max,
                y_p_max - leg_thick, y_p_max,
                -1.0, 2.75);

vase_body.setParams(4., 5., 3.8, 1.0);
vase_inner.setParams(4., 5., 3.3,
                     0., 0., 100.,
                     0.3);

CombiObjectOr desk (&plate);
desk.includeObject(&leg_1);
desk.includeObject(&leg_2);
desk.includeObject(&leg_3);
desk.includeObject(&leg_4);

CombiObjectAndNot vase (&vase_body);
vase.includeObject(&vase_inner);

CombiObjectOr object(&desk, &vase);

gridfile ="patches/standard.grid";
