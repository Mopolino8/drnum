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
#include "cylinderincartesianpatch.h"
#include "perfectgas.h"

CylinderInCartesianPatch::CylinderInCartesianPatch(CartesianPatch *patch)
{
  m_Patch  = patch;
  m_Temp   = 300.0;
  m_XOrg   = vec3_t(0.0, 0.0, 0.0);
  m_Rad    = 2.0;
  //m_Omega  = -1.04719775;
  m_Omega  = 17.3594;
}


void CylinderInCartesianPatch::operator ()()
{
  dim_t<5> dim;

  m_Patch->copyFieldToHost(0);

  real var1[5];
  for (int i_var = 0; i_var < 5; ++i_var) {
    var1[i_var] = 0;
  }
  // Clear residual field
  for (size_t i = 0; i < m_Patch->sizeI(); ++i) {
    for (size_t j = 0; j < m_Patch->sizeJ(); ++j) {
      for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
        m_Patch->setVar(dim, 2, i, j, k, var1);
      }
    }
  }
  vec3_t r, x_p;
  for (size_t i = 0; i < m_Patch->sizeI(); ++i) {
    for (size_t j = 0; j < m_Patch->sizeJ(); ++j) {
      for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
        x_p[0] = (real(i) + 0.5)*m_Patch->dx();
        x_p[1] = (real(j) + 0.5)*m_Patch->dy();
        x_p[2] = (real(k) + 0.5)*m_Patch->dz();
        x_p = m_Patch->getTransformInertial2This().transformReverse(x_p);

        r[0] = x_p[0] - m_XOrg[0];
        r[1] = x_p[1] - m_XOrg[1];
        r[2] = x_p[2] - m_XOrg[2];

        if ((r[0]*r[0] + r[1]*r[1]) < m_Rad*m_Rad) {
          // set to interior
          var1[0] = 1;
          // Abuse the residual field for storage
          m_Patch->setVar(dim, 2, i, j, k, var1);
          var1[0] = 0;
        }
        else {
          // clear data field for cells outside of cylinder
          m_Patch->setVar(dim, 0, i, j, k, var1);
        }
      }
    }
  }

  //
  real var2[5], d_var1[5], d_var2[5];
  vec3_t r_1, x_p1, r_2, x_p2;
  for (int i_var = 0; i_var < 5; ++i_var) {
    var1[i_var]   = 0;  // Residual field, pt 1
    var2[i_var]   = 0;  // Residual field, pt 2
    d_var1[i_var] = 0;  // Data field, pt 1
    d_var2[i_var] = 0;  // Data field, pt 2
  }
  for (size_t layer = 0; layer < 2; ++layer) {
    // iterate over I values, and compare i, i-1
    for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
      for (size_t i = 1; i < m_Patch->sizeI(); ++i) {
        for (size_t j = 0; j < m_Patch->sizeJ(); ++j) {
          m_Patch->getVar(dim, 2, i - 1, j, k, var1);
          m_Patch->getVar(dim, 2, i    , j, k, var2);

          // 0.5 to ensure real / int comparison works
          if ((var1[0] > (0.5 + layer) && var2[0] <  0.5 ) ||
              (var1[0] <  0.5          && var2[0] > (0.5 + layer)))
          {
            x_p1[0] = (real(i - 1) + 0.5)*m_Patch->dx();
            x_p1[1] = (real(j)     + 0.5)*m_Patch->dy();
            x_p1[2] = (real(k)     + 0.5)*m_Patch->dz();
            x_p1 = m_Patch->getTransformInertial2This().transformReverse(x_p1);
            r_1[0] = x_p1[0] - m_XOrg[0];
            r_1[1] = x_p1[1] - m_XOrg[1];
            r_1[2] = x_p1[2] - m_XOrg[2];

            x_p2[0] = (real(i) + 0.5)*m_Patch->dx();
            x_p2[1] = (real(j) + 0.5)*m_Patch->dy();
            x_p2[2] = (real(k) + 0.5)*m_Patch->dz();
            x_p2 = m_Patch->getTransformInertial2This().transformReverse(x_p2);
            r_2[0] = x_p2[0] - m_XOrg[0];
            r_2[1] = x_p2[1] - m_XOrg[1];
            r_2[2] = x_p2[2] - m_XOrg[2];

            m_Patch->getVar(dim, 0, i - 1, j, k, d_var1);
            m_Patch->getVar(dim, 0, i    , j, k, d_var2);

            if (var1[0] > (0.5 + layer)) {
              // pnt 1 is reference, interpolate to pnt 2
              real T, p, u, v, w;
              PerfectGas::conservativeToPrimitive(d_var1, p, T, u, v, w);

              real u_layer =  m_Omega*r_2[1];         // ==  m_Omega*norm(r_2)*r_2[1]/norm(r_2) ->  v_t*r*sin(th);
              real v_layer = -m_Omega*r_2[0];         // == -m_Omega*norm(r_2)*r_2[0]/norm(r_2) -> -v_t*r*cos(th);
              real weight  = scalarProduct(normalize(r_1), normalize(r_2));
              real dp_new  = weight*(p + m_Patch->dx()*(u*u + v*v + w*w)*d_var1[0]/norm(r_1));

              var2[0]  = -1;     // mark layer
              var2[1] += weight;
              var2[2] += dp_new;
              var2[3]  = u_layer;
              var2[4]  = v_layer;
              m_Patch->setVar(dim, 2, i, j, k, var2);
            } else {
              // pnt 2 is reference, extrapolate to pnt 1
              real T, p, u, v, w;
              PerfectGas::conservativeToPrimitive(d_var2, p, T, u, v, w);

              real u_layer =  m_Omega*r_1[1];         // == see above
              real v_layer = -m_Omega*r_1[0];         // == see above
              real weight  = scalarProduct(normalize(r_2), normalize(r_1));
              real dp_new  = weight*(p + m_Patch->dx()*(u*u + v*v + w*w)*d_var2[0]/norm(r_2));

              var1[0]  = -1;     // mark layer
              var1[1] += weight;
              var1[2] += dp_new;
              var1[3]  = u_layer;
              var1[4]  = v_layer;
              m_Patch->setVar(dim, 2, i - 1, j, k, var1);
            }
          } // end if loop
        } // end i loop
      } // end j loop
    } // end k loop

    // J values
    for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
      for (size_t i = 0; i < m_Patch->sizeI(); ++i) {
        for (size_t j = 1; j < m_Patch->sizeJ(); ++j) {
          m_Patch->getVar(dim, 2, i, j - 1, k, var1);
          m_Patch->getVar(dim, 2, i, j    , k, var2);

          if ((var1[0] > (0.5 + layer) && var2[0] <  0.5) ||
              (var1[0] <  0.5          && var2[0] > (0.5 + layer)))
          {
            x_p1[0] = (real(i)     + 0.5)*m_Patch->dx();
            x_p1[1] = (real(j - 1) + 0.5)*m_Patch->dy();
            x_p1[2] = (real(k)     + 0.5)*m_Patch->dz();
            x_p1 = m_Patch->getTransformInertial2This().transformReverse(x_p1);
            r_1[0] = x_p1[0] - m_XOrg[0];
            r_1[1] = x_p1[1] - m_XOrg[1];
            r_1[2] = x_p1[2] - m_XOrg[2];

            x_p2[0] = (real(i) + 0.5)*m_Patch->dx();
            x_p2[1] = (real(j) + 0.5)*m_Patch->dy();
            x_p2[2] = (real(k) + 0.5)*m_Patch->dz();
            x_p2 = m_Patch->getTransformInertial2This().transformReverse(x_p2);
            r_2[0] = x_p2[0] - m_XOrg[0];
            r_2[1] = x_p2[1] - m_XOrg[1];
            r_2[2] = x_p2[2] - m_XOrg[2];

            m_Patch->getVar(dim, 0, i, j - 1, k, d_var1);
            m_Patch->getVar(dim, 0, i, j    , k, d_var2);

            if (var1[0] > (0.5 + layer) && var2[0] < (0.5 + layer)) {
              // pnt 1 is reference, extrapolate to pnt 2
              real T, p, u, v, w;
              PerfectGas::conservativeToPrimitive(d_var1, p, T, u, v, w);

              real u_layer =  m_Omega*r_2[1];         // ==  m_Omega*norm(r_2)*r_2[1]/norm(r_2) ->  v_t*r*sin(th);
              real v_layer = -m_Omega*r_2[0];         // == -m_Omega*norm(r_2)*r_2[0]/norm(r_2) -> -v_t*r*cos(th);
              real weight  = scalarProduct(normalize(r_1), normalize(r_2));
              real dp_new  = weight*(p + m_Patch->dy()*(u*u + v*v + w*w)*d_var1[0]/norm(r_1));

              var2[0]  = -1;     // mark layer
              var2[1] += weight;
              var2[2] += dp_new;
              var2[3]  = u_layer;
              var2[4]  = v_layer;
              m_Patch->setVar(dim, 2, i, j, k, var2);
            }
            else if (var1[0] < (0.5 + layer) && var2[0] > (0.5 + layer)) {
              // pnt 2 is reference, extrapolate to pnt 1
              real T, p, u, v, w;
              PerfectGas::conservativeToPrimitive(d_var2, p, T, u, v, w);

              real u_layer =  m_Omega*r_1[1];         // == see above
              real v_layer = -m_Omega*r_1[0];         // == see above
              real weight = scalarProduct(normalize(r_2), normalize(r_1));
              real dp_new = weight*(p + m_Patch->dy()*(u*u + v*v + w*w)*d_var2[0]/norm(r_2));

              var1[0] = -1;     // mark layer
              var1[1] += weight;
              var1[2] += dp_new;
              var1[3] = u_layer;
              var1[4] = v_layer;
              m_Patch->setVar(dim, 2, i, j - 1, k, var1);
            }
          } // end if loop
        } // end i loop
      } // end j loop
    } // end k loop

    // Set the marked layers to the correct type number. 2 for layer 1, 3 for layer 2.
    // Then update the data field
    for (size_t i = 0; i < m_Patch->sizeI(); ++i) {
      for (size_t j = 1; j < m_Patch->sizeJ(); ++j) {
        for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
          m_Patch->getVar(dim, 2, i, j, k, var1);
          if (var1[0] < -0.5) {
            var1[0] = 2 + layer;
            m_Patch->setVar(dim, 2, i, j, k, var1);

            m_Patch->getVar(dim, 0, i, j, k, d_var1);
            real u, v, p, weight;
            weight = var1[1];
            p      = var1[2];
            u      = var1[3];
            v      = var1[4];
            PerfectGas::primitiveToConservative(p/weight, m_Temp, u, v, 0, d_var1);
            m_Patch->setVar(dim, 0, i, j, k, d_var1);
          }
        } // end i loop
      } // end j loop
    } // end k loop

  } // end layer loop


  // Set the outside field to zero.
  for (size_t i = 0; i < m_Patch->sizeI(); ++i) {
    for (size_t j = 0; j < m_Patch->sizeJ(); ++j) {
      for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
        m_Patch->getVar(dim, 2, i, j, k, var1);
        x_p2[0] = (real(i) + 0.5)*m_Patch->dx();
        x_p2[1] = (real(j) + 0.5)*m_Patch->dy();
        x_p2[2] = (real(k) + 0.5)*m_Patch->dz();
        x_p2 = m_Patch->getTransformInertial2This().transformReverse(x_p2);
        if (var1[0] < 0.5) {
          PerfectGas::primitiveToConservative(1e5, m_Temp, 0, 0, 0, d_var1);
          m_Patch->setVar(dim, 0, i, j, k, d_var1);
        }
          /*
        if (var1[0] < 1.5 && var1[0] > 0.5) {
          cout << "data1 " << x_p2[0]
               << " " << x_p2[1]
               << " " << x_p2[2]
               << endl;
        }
        if (var1[0] < 2.5 && var1[0] > 1.5) {
          cout << "data2 " << x_p2[0]
               << " " << x_p2[1]
               << " " << x_p2[2]
               << endl;
        }
        if (var1[0] > 2.5) {
          cout << "data3 " << x_p2[0]
               << " " << x_p2[1]
               << " " << x_p2[2]
               << endl;
        }
        */
      } // end i loop
    } // end j loop
  } // end k loop

  // Clear residual field
  for (int i_var = 0; i_var < 5; ++i_var) {
    var1[i_var] = 0;
  }
  for (size_t i = 0; i < m_Patch->sizeI(); ++i) {
    for (size_t j = 0; j < m_Patch->sizeJ(); ++j) {
      for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
        m_Patch->setVar(dim, 2, i, j, k, var1);
      }
    }
  }

  m_Patch->copyFieldToDevice(0);
}
