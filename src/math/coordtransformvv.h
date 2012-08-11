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
