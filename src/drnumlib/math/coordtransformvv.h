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
#ifndef COORDTRANSFORMVV_H
#define COORDTRANSFORMVV_H

#include "blockcfd.h"
#include "coordtransform.h"

class CoordTransformVV
{
  CoordTransform m_transform;
  CoordTransform m_transform_inverse;

private:

  void computeInverse();

public:

  CoordTransformVV();

  void setMatrix(const mat3_t& A);

  /**
    * Set up turn matrix to transform base vectors.
    * (1,0,0) into normal(base_i)
    * (0,1,0) into normal(base_j - dotprod(base_i, base_j) * base_i)
    * (0,0,1) into dotprod of above two
    *
    * NOTES:
    *  *base vector af j-coords will be modified to get dotprod(i_base, j_base) == 0
    *  *base vector of k-coordinates is computed as crossprod(i_base, j_base)
    * @param base_i the base vector of i-coordinates in non transformed system
    * @param base_j the base vector of j-coordinates in non transformed system
    * @return abs(crossprod(norm(i_base), norm(j_base))) . Approaches zero for lin. dependency.
    */
  real setMatrixFromBaseIJ(vec3_t base_i, vec3_t base_j);

  void setVector(const vec3_t& b);
  void setAll(const mat3_t& A, const vec3_t& b);
  void setAll(const CoordTransform& transform);
  void setTransFromTo(const CoordTransformVV& c_from, const CoordTransformVV& c_to);
  mat3_t getMatrix();
  vec3_t getVector();
  mat3_t getInvMatrix();
  vec3_t getInvVector();
  CoordTransform extractForward();
  CoordTransform extractReverse();

  /**
    * Scale shift vector.
    * @param scfactor the scaling factor.
    */
  void scaleVector(real scfactor);

  CoordTransformVV concatenate(const CoordTransformVV& c) const;  ///< ret = this * c;
  void concatenate_reflexive(const CoordTransformVV& c);          ///< this = this * c;

  void operator=(const CoordTransformVV& c);  /// @todo needed? Try standard copy.
  void copy(const CoordTransformVV& c) {operator=(c);}

  void transform(real x, real y, real z, real& xt, real& yt, real& zt) const;
  void transform(real& x, real& y, real& z) const { transform(x, y, z, x, y, z); }
  vec3_t transform(const vec3_t& xyz) const;

  void transformReverse(const real& x, const real& y, const real& z,
                        real& xt, real& yt, real& zt) const;
  void transformReverse(real& x, real& y, real& z) const {
    transformReverse(x, y, z, x, y, z);}
  vec3_t transformReverse(const vec3_t& xyz) const;

  void transfree(const real& u, const real& v, const real& w,
                 real& ut, real& vt, real& wt) const;
  void transfree(real& u, real& v, real& w) const {
    transfree(u, v, w, u, v, w);}
  vec3_t transfree(const vec3_t& uvw) const;

  void transfreeReverse(const real& u, const real& v, const real& w,
                        real& ut, real& vt, real& wt) const;
  void transfreeReverse(real& u, real& v, real& w) const {
    transfreeReverse(u, v, w, u, v, w);}
  vec3_t transfreeReverse(const vec3_t& uvw) const;

};


inline void CoordTransformVV::operator=(const CoordTransformVV& c)
{
  m_transform = c.m_transform;
  m_transform_inverse = c.m_transform_inverse;
}


inline void CoordTransformVV::transform(real x, real y, real z, real& xt, real& yt, real& zt) const
{
  m_transform.transform(x, y, z, xt, yt, zt);
}

inline vec3_t CoordTransformVV::transform(const vec3_t& xyz) const
{
  vec3_t xyzt;
  m_transform.transform(xyz[0], xyz[1], xyz[2], xyzt[0], xyzt[1], xyzt[2]);
  return xyzt;
}

inline void CoordTransformVV::transformReverse(const real& x, const real& y, const real& z,
                                               real& xt, real& yt, real& zt) const
{
  m_transform_inverse.transform(x, y, z, xt, yt, zt);
}

inline vec3_t CoordTransformVV::transformReverse(const vec3_t& xyz) const
{
  vec3_t xyzt;
  m_transform_inverse.transform(xyz[0], xyz[1], xyz[2], xyzt[0], xyzt[1], xyzt[2]);
  return xyzt;
}

inline void CoordTransformVV::transfree(const real& u, const real& v, const real& w,
                                        real& ut, real& vt, real& wt) const
{
  m_transform.transfree(u, v, w, ut, vt, wt);
}

inline vec3_t CoordTransformVV::transfree(const vec3_t& uvw) const
{
  vec3_t uvwt;
  m_transform.transfree(uvw[0], uvw[1], uvw[2], uvwt[0], uvwt[1], uvwt[2]);
  return uvwt;
}

inline void CoordTransformVV::transfreeReverse(const real& u, const real& v, const real& w,
                                               real& ut, real& vt, real& wt) const
{
  m_transform_inverse.transfree(u, v, w, ut, vt, wt);
}

inline vec3_t CoordTransformVV::transfreeReverse(const vec3_t& uvw) const
{
  vec3_t uvwt;
  m_transform_inverse.transfree(uvw[0], uvw[1], uvw[2], uvwt[0], uvwt[1], uvwt[2]);
  return uvwt;
}

inline CoordTransformVV CoordTransformVV::concatenate(const CoordTransformVV &c) const
{
  CoordTransformVV t;
  t.m_transform = m_transform.concatenate(c.m_transform);
  t.computeInverse();
  return t;
}

inline void CoordTransformVV::concatenate_reflexive(const CoordTransformVV &c)
{
  // reflexive version of concatenate
  m_transform.concatenate_reflexive(c.m_transform);
  computeInverse();
}



#endif // COORDTRANSFORMVV_H
