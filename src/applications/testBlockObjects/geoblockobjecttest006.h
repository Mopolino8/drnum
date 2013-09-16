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
// + DrNUM is distributed in the hope that it will be useful,             +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with DrNUM. If not, see <http://www.gnu.org/licenses/>.        +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "blockcfd.h"

// Includes to build a geometric test case

// rocket example needs grid rocket.grid

real x_base = 0.;
real x_stage_1 = 1.;
real x_stage_2 = 2.8;
real x_nose_cone = 3.5;
real x_nose = 4.0;
real x_boost_nose = 1.5;

real y_c = 0.;
real y_left = -0.4;
real y_right = 0.4;

real z_all = 0.;

real r_thick = 0.25;
real r_thin = 0.15;
real r_boost = 0.2;
real r_min = 0.03;

real intercone = 0.2;

real save = 0.01;

// stage 1;
CylinderObject stage_1_cyl;
stage_1_cyl.setParams(x_base, y_c, z_all,
                      (x_stage_1 - x_base), 0., 0.,
                      r_thick);

// stage 2:
CylinderObject stage_2_cyl;
stage_2_cyl.setParams((x_stage_1 - save), y_c, z_all,
                      (x_stage_2 - x_stage_1 + save), 0., 0.,
                      r_thin);

// stage 3:
CylinderObject stage_3_cyl;
stage_3_cyl.setParams(x_stage_2, y_c, z_all,
                      (x_nose_cone - x_stage_2), 0., 0.,
                      r_thick);

// interstage 1-2:
ConeObject interstage_1_2;
interstage_1_2.setParams(x_stage_1, y_c, z_all,
                         (2*intercone), 0., 0.,
                         r_thick, r_thin);

// interstage 2-3:
ConeObject interstage_2_3;
interstage_2_3.setParams(x_stage_2, y_c, z_all,
                         - intercone, 0., 0.,
                         r_thick, r_thin);

// nose cone
ConeObject nose_cone;
nose_cone.setParams(x_nose_cone, y_c, z_all,
                    x_nose - x_nose_cone, 0., 0.,
                    r_thick, r_min);

// left booster
//.. body
CylinderObject lb_cyl;
lb_cyl.setParams(x_base, y_left, z_all,
                 (x_stage_1 - x_base), 0., 0.,
                 r_boost);
//.. cone
ConeObject lb_cone;
lb_cone.setParams(x_stage_1, y_left, z_all,
                  (x_boost_nose - x_stage_1),  0., 0.,
                  r_boost, r_min);

// right booster
//.. body
CylinderObject rb_cyl;
rb_cyl.setParams(x_base, y_right, z_all,
                 (x_stage_1 - x_base), 0., 0.,
                 r_boost);
//.. cone
ConeObject rb_cone;
rb_cone.setParams(x_stage_1, y_right, z_all,
                  (x_boost_nose - x_stage_1),  0., 0.,
                  r_boost, r_min);

// all together
CombiObjectOr object (&stage_1_cyl);
object.includeObject (&stage_2_cyl);
object.includeObject (&stage_3_cyl);
object.includeObject (&interstage_1_2);
object.includeObject (&interstage_2_3);
object.includeObject (&nose_cone);
object.includeObject (&lb_cyl);
object.includeObject (&lb_cone);
object.includeObject (&rb_cyl);
object.includeObject (&rb_cone);

gridfile = "patches/rocket.grid";
