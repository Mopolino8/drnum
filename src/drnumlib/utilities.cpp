// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
// +                                                                      +
// + enGrid is free software: you can redistribute it and/or modify       +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 
#include "utilities.h"

#include <cmath>

#include <QString>
#include <QStringList>
#include <QDir>

//#include "error.h"
#include "math/mathvector.h"
#include "math/smallsquarematrix.h"
#include "math/linsolve.h"


QString addSuffix(QString filename, QString suffix, bool remove_old_suffix) {
  QFileInfo fileinfo(filename);
  if (remove_old_suffix) {
    return fileinfo.completeBaseName() + QObject::tr(".") + suffix;
  } else {
    if (fileinfo.suffix().toLower() == suffix) {
      return fileinfo.absoluteFilePath();
    } else {
      return fileinfo.absoluteFilePath() + QObject::tr(".") + suffix;
    }
  }
}

mat3_t rotationMatrixX(real a_rad)
{
  mat3_t Rx;
  Rx[0][0] = 1; Rx[0][1] = 0;           Rx[0][2] = 0;
  Rx[1][0] = 0; Rx[1][1] = cos(a_rad);  Rx[1][2] = sin(a_rad);
  Rx[2][0] = 0; Rx[2][1] = -sin(a_rad); Rx[2][2] = cos(a_rad);
  return Rx;
}

mat3_t rotationMatrixY(real a_rad)
{
  mat3_t Ry;
  Ry[0][0] = cos(a_rad); Ry[0][1] = 0; Ry[0][2] = -sin(a_rad);
  Ry[1][0] = 0;          Ry[1][1] = 1; Ry[1][2] = 0;
  Ry[2][0] = sin(a_rad); Ry[2][1] = 0; Ry[2][2] = cos(a_rad);
  return Ry;
}

mat3_t rotationMatrixZ(real a_rad)
{
  mat3_t Rz;
  Rz[0][0] = cos(a_rad);  Rz[0][1] = sin(a_rad); Rz[0][2] = 0;
  Rz[1][0] = -sin(a_rad); Rz[1][1] = cos(a_rad); Rz[1][2] = 0;
  Rz[2][0] = 0;           Rz[2][1] = 0;          Rz[2][2] = 1;
  return Rz;
}

mat3_t rotateRelativeZXZ(real angle_1_rad, real angle_2_rad, real angle_3_rad)
{
  return rotationMatrixZ(angle_3_rad)*rotationMatrixX(angle_2_rad)*rotationMatrixZ(angle_1_rad);
}

mat3_t rotateAbsoluteZXY(real angle_1_rad, real angle_2_rad, real angle_3_rad)
{
  return rotationMatrixZ(angle_1_rad)*rotationMatrixX(angle_2_rad)*rotationMatrixY(angle_3_rad);
}

real getGamma(vec3_t V) {
  return V[0] == 0.0 && V[1] == 0.0 && V[2] == 0.0 ? 0.0 : atan2(sqrt(V[0]*V[0] + V[1]*V[1]), V[2]);
}

real getPhi(vec3_t V) {
  return V[0] == 0.0 && V[1] == 0.0 ? 0.0 : atan2(V[1], V[0]);
}

bool checkVector(vec3_t V)
{
  for (int i = 0; i < 3; i++) {
    if (isnan(V[i])) {
      return false;
    }
    if (isinf(V[i])) {
      return false;
    }
  }
  return true;
}

bool checkVector(vec2_t V)
{
  for (int i = 0; i < 2; i++) {
    if (isnan(V[i])) {
      return false;
    }
    if (isinf(V[i])) {
      return false;
    }
  }
  return true;
}

int polySolveCubic(real a, real b, real c, real * x0, real * x1, real * x2)
{
  real q = (a * a - 3 * b);
  real r = (2 * a * a * a - 9 * a * b + 27 * c);
  
  real Q = q / 9;
  real R = r / 54;
  
  real Q3 = Q * Q * Q;
  real R2 = R * R;
  
  real CR2 = 729 * r * r;
  real CQ3 = 2916 * q * q * q;
  
  if (R == 0 && Q == 0)
  {
    *x0 = - a / 3 ;
    *x1 = - a / 3 ;
    *x2 = - a / 3 ;
    return 3 ;
  }
  else if (CR2 == CQ3) 
  {
    // this test is actually R2 == Q3, written in a form suitable
    // for exact computation with integers
    
    // Due to finite precision some real roots may be missed, and
    // considered to be a pair of complex roots z = x +/- epsilon i
    // close to the real axis.
    
    real sqrtQ = sqrt (Q);
    
    if (R > 0)
    {
      *x0 = -2 * sqrtQ  - a / 3;
      *x1 = sqrtQ - a / 3;
      *x2 = sqrtQ - a / 3;
    } else {
      *x0 = - sqrtQ  - a / 3;
      *x1 = - sqrtQ - a / 3;
      *x2 = 2 * sqrtQ - a / 3;
    }
    return 3;
  } else if (CR2 < CQ3) { // equivalent to R2 < Q3
    real sqrtQ = sqrt (Q);
    real sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
    real theta = acos (R / sqrtQ3);
    real norm = -2 * sqrtQ;
    *x0 = norm * cos (theta / 3) - a / 3;
    *x1 = norm * cos ((theta + 2.0 * M_PI) / 3) - a / 3;
    *x2 = norm * cos ((theta - 2.0 * M_PI) / 3) - a / 3;
    
    // Sort *x0, *x1, *x2 into increasing order
    
    if (*x0 > *x1) {
      swap(*x0, *x1);
    }
    if (*x1 > *x2)
    {
      swap(*x1, *x2);
      if (*x0 > *x1) {
        swap(*x0, *x1);
      }
    }
    
    return 3;
  } else {
    real sgnR = (R >= 0 ? 1 : -1);
    real A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), 1.0/3.0);
    real B = Q / A ;
    *x0 = A + B - a / 3;
    return 1;
  }
}

// a x^2 + b x + c = 0
int poly_solve_quadratic(real a, real b, real c, real * x0, real * x1)
{
  if (a == 0) {
    if (b == 0) {
      return(0);
    } else {
      *x0 = -c / b;
      return(1);
    }
  } else {
    real delta = pow(b, 2) - 4 * a * c;
    if (delta < 0) {
      return(0);
    } else {
      *x0 = (-b + sqrt(delta)) / (2 * a);
      *x1 = (-b - sqrt(delta)) / (2 * a);
      if (*x1 < *x0) {
        real tmp = *x0;
        *x0 = *x1;
        *x1 = tmp;
      }
    }
  }
  return(2);
}
