#ifndef COORDTRANSFORM_H
#define COORDTRANSFORM_H

#include "blockcfd.h"

class CoordTransform
{
  real axx, axy, axz, bx;
  real ayx, ayy, ayz, by;
  real azx, azy, azz, bz;


public:

  CoordTransform();

  void setMatrix(mat3_t A);

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

  void setVector(vec3_t b);
  void setAll(mat3_t A, vec3_t b);
  void setTransFromTo(const CoordTransform& c_from, const CoordTransform& c_to);
  mat3_t getMatrix();
  vec3_t getVector();
  CoordTransform inverse() const;

  /**
    * Scale shift vector.
    * @param scfactor the scaling factor.
    */
  void scaleVector(real scfactor);

  CoordTransform concatenate(const CoordTransform& c) const;  ///< ret = this * c;
  /// @todo cant overload
  void concatenate_reflexive(const CoordTransform& c);        ///< this = this * c;

  //void operator=(const CoordTransform& c);
  void copy(const CoordTransform& c) {operator=(c);}

  void transform(real x, real y, real z, real& xt, real& yt, real& zt) const;
  void transform(real& x, real& y, real& z) const { transform(x, y, z, x, y, z); }
  vec3_t transform(const vec3_t& xyz) const;

  void transfree(const real& u, const real& v, const real& w,
                 real& ut, real& vt, real& wt) const;
  void transfree(real& u, real& v, real& w) const {
    transfree(u, v, w, u, v, w);}
  vec3_t transfree(const vec3_t& uvw) const;

};


//void CoordTransform::operator=(const CoordTransform& c)
//{
//  axx = c.axx; axy = c.axy; axz = c.axz; bx = c.bx;
//  ayx = c.ayx; ayy = c.ayy; ayz = c.ayz; by = c.by;
//  azx = c.azx; azy = c.azy; azz = c.azz; bz = c.bz;
//}


inline void CoordTransform::transform(real x, real y, real z, real& xt, real& yt, real& zt) const
{
  xt = axx*x + axy*y + axz*z + bx;
  yt = ayx*x + ayy*y + ayz*z + by;
  zt = azx*x + azy*y + azz*z + bz;
  countFlops(18);
}

inline vec3_t CoordTransform::transform(const vec3_t& xyz) const
{
  vec3_t xyzt;
  transform(xyz[0], xyz[1], xyz[2], xyzt[0], xyzt[1], xyzt[2]);
  return xyzt;
}

inline void CoordTransform::transfree(const real& u, const real& v, const real& w,
                                      real& ut, real& vt, real& wt) const
{
  ut = axx*u + axy*v + axz*w;
  vt = ayx*u + ayy*v + ayz*w;
  wt = azx*u + azy*v + azz*w;
  countFlops(15);
}

inline vec3_t CoordTransform::transfree(const vec3_t& uvw) const
{
  vec3_t uvwt;
  transfree(uvw[0], uvw[1], uvw[2], uvwt[0], uvwt[1], uvwt[2]);
  return uvwt;
}

inline CoordTransform CoordTransform::concatenate(const CoordTransform &c) const
{
//inline CoordTransform CoordTransform::concatenate(const CoordTransform &c) const
//{
  /** @todo Check carefully. Compared to transformation.h
   * I've eliminated all multiplications involving 0.
   * There is no need to store fourth row, as it allways reads (0 0 0 1) */

  CoordTransform t;

  t.axx = axx*c.axx + axy*c.ayx + axz*c.azx;
  t.axy = axx*c.axy + axy*c.ayy + axz*c.azy;
  t.axz = axx*c.axz + axy*c.ayz + axz*c.azz;
  t.bx = axx*c.bx + axy*c.by + axz*c.bz + bx;
  countFlops(21);

  t.ayx = ayx*c.axx + ayy*c.ayx + ayz*c.azx;
  t.ayy = ayx*c.axy + ayy*c.ayy + ayz*c.azy;
  t.ayz = ayx*c.axz + ayy*c.ayz + ayz*c.azz;
  t.by = ayx*c.bx + ayy*c.by + ayz*c.bz + by;
  countFlops(21);

  t.azx = azx*c.axx + azy*c.ayx + azz*c.azx;
  t.azy = azx*c.axy + azy*c.ayy + azz*c.azy;
  t.azz = azx*c.axz + azy*c.ayz + azz*c.azz;
  t.bz = azx*c.bx + azy*c.by + azz*c.bz + bz;
  countFlops(21);

  return t;
}

inline void CoordTransform::concatenate_reflexive(const CoordTransform &c)
{
  // reflexive version of concatenate
  copy(concatenate(c));
}

#endif // COORDTRANSFORM_H
