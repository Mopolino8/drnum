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
#include "geometrytools.h"
#include "containertricks.h"
//#include "utilities.h"
#include "math/linsolve.h"

#include <vtkCellType.h>
#include <cmath>

namespace GeometryTools
{

real rad2deg(real rad)
{
  return rad/M_PI*180;
}

real deg2rad(real deg)
{
  return deg/180*M_PI;
}

void rotate(vec3_t g1, vec3_t g2,vec3_t g3, vec3_t &b, real theta)
{
  //cout << b << endl;
  mat3_t g_t;
  g_t[0] = g1;
  g_t[1] = g2;
  g_t[2] = g3;
  mat3_t g = g_t.transp();
  //cout << g << endl;
  vec3_t bb = g_t*b;
  //cout << bb << endl;
  mat3_t rot = mat3_t::identity();
  rot[0][0] = cos(theta);
  rot[0][1] = -sin(theta);
  rot[1][0] = sin(theta);
  rot[1][1] = cos(theta);
  //cout << rot << endl;
  vec3_t bbb = rot*bb;
  //cout << bbb << endl;
  b = g*bbb;
  //cout << b << endl;
}

vec3_t rotate(vec3_t v, vec3_t axis, real theta)
{
  axis.normalise();
  
  // transposed base of rotate system
  mat3_t g_t;

  // compute projection of v in axis direction
  vec3_t v_axis = (axis*v)*axis;

  // compute first orthogonal vector (first base vector)
  g_t[0] = v-v_axis;
  
  //In case of points on the rotation axis, do nothing
  if(g_t[0].abs()==0) return v;
  
  g_t[0].normalise();

  // second base vector is the normalised axis
  g_t[1] = axis;

  // compute second orthogonal vector (third base vector)
  g_t[2] = g_t[0].cross(g_t[1]);

  // base of rotate system
  mat3_t g = g_t.transp();

  // matrix for rotation around g_t[1];
  mat3_t rot = mat3_t::identity();
  rot[0][0] =  cos(theta);
  rot[0][2] =  sin(theta);
  rot[2][0] = -sin(theta);
  rot[2][2] =  cos(theta);

  // transfer v to rotate system
  vec3_t v_r = g_t*v;

  // rotate the vector and transfer it back
  v_r = rot*v_r;
  v = g*v_r;

  return v;
}

vec3_t orthogonalVector(vec3_t v)
{
  // get absolute values
  real xx = v[0] < 0.0 ? -v[0] : v[0];
  real yy = v[1] < 0.0 ? -v[1] : v[1];
  real zz = v[2] < 0.0 ? -v[2] : v[2];
  // switch both biggest values and set the other one to zero
  vec3_t u;
  if (xx < yy) {
    u = xx < zz ? vec3_t(0,v[2],-v[1]) : vec3_t(v[1],-v[0],0);
  } else {
    u = yy < zz ? vec3_t(-v[2],0,v[0]) : vec3_t(v[1],-v[0],0);
  }
  u.normalise();
  return u;
}


real intersection(vec3_t x_straight, vec3_t v_straight, vec3_t x_plane, vec3_t u_plane, vec3_t v_plane)
{
  vec3_t n = u_plane.cross(v_plane);
  return intersection(x_straight,v_straight,x_plane,n);
}

real intersection(vec3_t x_straight, vec3_t v_straight, vec3_t x_plane, vec3_t n_plane)
{
  real k = (x_plane*n_plane - x_straight*n_plane)/(v_straight*n_plane);
  return k;
}

bool intersection (real &k1, real &k2, vec2_t r1, vec2_t u1, vec2_t r2, vec2_t u2)
{
  real ave_length = .5*(sqrt(u1[0]*u1[0]+u1[1]*u1[1]) + sqrt(u2[0]*u2[0]+u2[1]*u2[1]));
  real DET = (u1[0]*u2[1]-u1[1]*u2[0]);
  if (fabs(DET) > 1e-6*ave_length) {
    k1 = -(u2[0]*r2[1]-u2[0]*r1[1]-r2[0]*u2[1]+r1[0]*u2[1])/DET;
    k2 = -(-u1[1]*r2[0]+u1[0]*r2[1]-u1[0]*r1[1]+u1[1]*r1[0])/DET;
    return true;
  } else {
    return false;
  }
}

void sliceTriangle(const vector<vec3_t> &Tin, vec3_t x, vec3_t n, vector<vector<vec3_t> > &Tout)
{
  vec3_t a = Tin[0];
  vec3_t b = Tin[1];
  vec3_t c = Tin[2];
  real kab = intersection(a,b-a,x,n);
  real kbc = intersection(b,c-b,x,n);
  real kca = intersection(c,a-c,x,n);
  bool ab_cut = ((kab >= 0) && (kab <= 1));
  bool bc_cut = ((kbc >= 0) && (kbc <= 1));
  bool ca_cut = ((kca >= 0) && (kca <= 1));
  if (ab_cut && bc_cut && ca_cut) {
    //cerr << "invalid triangle (SliceTriangle) A" << endl;
    //exit(EXIT_FAILURE);
    if      ((kab <= kbc) && (kab <= kca)) ab_cut = false;
    else if ((kbc <= kab) && (kbc <= kca)) bc_cut = false;
    else                                   ca_cut = false;
  }
  if (ab_cut && bc_cut) {
    vec3_t ab = a + kab*(b-a);
    vec3_t bc = b + kbc*(c-b);
    Tout.resize(3,vector<vec3_t>(3));
    clinit(Tout[0]) = a,ab,bc;
    clinit(Tout[1]) = ab,b,bc;
    clinit(Tout[2]) = bc,c,a;
  } else if (bc_cut && ca_cut) {
    vec3_t bc = b + kbc*(c-b);
    vec3_t ca = c + kca*(a-c);
    Tout.resize(3,vector<vec3_t>(3));
    clinit(Tout[0]) = a,bc,ca;
    clinit(Tout[1]) = a,b,bc;
    clinit(Tout[2]) = bc,c,ca;
  } else if (ca_cut && ab_cut) {
    vec3_t ca = c + kca*(a-c);
    vec3_t ab = a + kab*(b-a);
    Tout.resize(3,vector<vec3_t>(3));
    clinit(Tout[0]) = a,ab,ca;
    clinit(Tout[1]) = ab,b,ca;
    clinit(Tout[2]) = b,c,ca;
  } else {
    Tout.resize(1,vector<vec3_t>(3));
    clinit(Tout[0]) = a,b,c;
  }
}

real tetraVol(const vec3_t& x0, const vec3_t& x1, const vec3_t& x2, const vec3_t& x3, bool neg)
{
  static real f16 = 1.0/6.0;
  vec3_t v1(x1[0]-x0[0], x1[1]-x0[1], x1[2]-x0[2]);
  vec3_t v2(x2[0]-x0[0], x2[1]-x0[1], x2[2]-x0[2]);
  vec3_t v3(x3[0]-x0[0], x3[1]-x0[1], x3[2]-x0[2]);
  real V = v1[0]*(v2[1]*v3[2]-v2[2]*v3[1]) + v1[1]*(v2[2]*v3[0]-v2[0]*v3[2]) + v1[2]*(v2[0]*v3[1]-v2[1]*v3[0]);
  V *= f16;
  if (!neg && (V < 0)) {
    V = -1e99;
  }
  return V; //fabs(V);
}

real pyraVol(vec3_t x0, vec3_t x1, vec3_t x2, vec3_t x3, vec3_t x4, bool neg)
{
  real V1 = tetraVol(x0, x1, x3, x4, neg) + tetraVol(x1, x2, x3, x4, neg);
  real V2 = tetraVol(x0, x1, x2, x4, neg) + tetraVol(x2, x3, x0, x4, neg);
  return min(V1,V2);
  /*
  real V = 0;
  vec3_t m0 = .25*(x0+x1+x2+x3);
  V += tetraVol(x0, x1, m0, x4, neg);
  V += tetraVol(x1, x2, m0, x4, neg);
  V += tetraVol(x2, x3, m0, x4, neg);
  V += tetraVol(x3, x0, m0, x4, neg);
  return V;
  */
}

real prismVol(vec3_t x0, vec3_t x1, vec3_t x2, vec3_t x3, vec3_t x4, vec3_t x5, bool neg)
{
  real V = 0;
  vec3_t p  = 1.0/6.0*(x0+x1+x2+x3+x4+x5);
  V += tetraVol(x0, x2, x1, p, neg);
  V += tetraVol(x3, x4, x5, p, neg);
  V += pyraVol (x0, x1, x4, x3, p, neg);
  V += pyraVol (x1, x2, x5, x4, p, neg);
  V += pyraVol (x0, x3, x5, x2, p, neg);
  return V;
}

real hexaVol(vec3_t x0, vec3_t x1, vec3_t x2, vec3_t x3, vec3_t x4, vec3_t x5, vec3_t x6, vec3_t x7, bool neg)
{
  real V = 0;
  vec3_t p = 1.0/8.0*(x0+x1+x2+x3+x4+x5+x6+x7);
  V += pyraVol(x0, x1, x3, x2, p, neg);
  V += pyraVol(x0, x4, x5, x1, p, neg);
  V += pyraVol(x4, x6, x7, x5, p, neg);
  V += pyraVol(x2, x3, x7, x6, p, neg);
  V += pyraVol(x1, x5, x7, x3, p, neg);
  V += pyraVol(x0, x2, x6, x4, p, neg);
  return V;
}

real triArea(vec3_t x0, vec3_t x1, vec3_t x2)
{
  vec3_t a = x1-x0;
  vec3_t b = x2-x0;
  real A = 0.5*((a.cross(b)).abs());
  return A;
}

real quadArea(vec3_t x0, vec3_t x1, vec3_t x2, vec3_t x3)
{
  real A = 0;
  vec3_t p = .25*(x0+x1+x2+x3);
  A += triArea(x0,x1,p);
  A += triArea(x1,x2,p);
  A += triArea(x2,x3,p);
  A += triArea(x3,x0,p);
  return A;
}

vec3_t triNormal(vec3_t x0, vec3_t x1, vec3_t x2)
{
  vec3_t a = x1-x0;
  vec3_t b = x2-x0;
  vec3_t n = 0.5*(a.cross(b));
  return n;
}

vec3_t quadNormal(vec3_t x0, vec3_t x1, vec3_t x2, vec3_t x3)
{
  vec3_t n;
  clinit(n) = 0,0,0;
  vec3_t p = .25*(x0+x1+x2+x3);
  n += triNormal(x0,x1,p);
  n += triNormal(x1,x2,p);
  n += triNormal(x2,x3,p);
  n += triNormal(x3,x0,p);
  return n;
}

real angle(const vec3_t & u, const vec3_t & v)
{
  // return the angle w.r.t. another 3-vector
  real ptot2 = u.abs2()*v.abs2();
  if(ptot2 <= 0) {
      return 0.0;
  } else {
    real arg = (u*v)/sqrt(ptot2);
    if(arg >  1.0) arg =  1.0;
    if(arg < -1.0) arg = -1.0;
    return acos(arg);
  }
}

bool isInsideTriangle(vec2_t r, real tol)
{
  if (r[0] < 0-tol || 1+tol < r[0] || r[1] < 0-tol || 1+tol < r[1] || r[0]+r[1] > 1+tol) {
    return false;
  }
  return true;
}

bool intersectEdgeAndTriangle(const vec3_t& a, const vec3_t& b, const vec3_t& c,
                              const vec3_t& x1, const vec3_t& x2, vec3_t& xi, vec3_t& ri, real tol)
{
  // triangle base
  vec3_t g1 = b - a;
  vec3_t g2 = c - a;
  vec3_t g3 = g1.cross(g2);
  g3.normalise();

  // direction of the edge
  vec3_t v = x2 - x1;

  // parallel?
  if (fabs(g3*v) < 1e-6) {
    return false;
  }

  // compute intersection between straight and triangular plane
  real k = intersection(x1, v, a, g3);
  xi = x1 + k*v;

  // transform xi to triangular base
  mat3_t G;
  G.column(0, g1);
  G.column(1, g2);
  G.column(2, g3);
  
  mat3_t GI = G.inverse();
  ri = xi - a;
  ri = GI*ri;
  /*
  {
    vec3_t b = xi-a;
    linsolve(G, b, ri);
  }
  */

  // intersection outside of edge range?
  if (k < 0 - tol) {
    return false;
  }
  if (k > 1 + tol) {
    return false;
  }
  
  // intersection outside of triangle?
  if (!isInsideTriangle(vec2_t(ri[0],ri[1]),tol)) {
    return false;
  }
  return true;
}

void computeCircumscribedCircle(vec3_t a, vec3_t b, vec3_t c, vec3_t &x, real &radius)
{
  real la = (b-c).abs();
  real lb = (a-c).abs();
  real lc = (a-b).abs();
  real bca = sqr(la)*(sqr(lb) + sqr(lc) - sqr(la));
  real bcb = sqr(lb)*(sqr(la) + sqr(lc) - sqr(lb));
  real bcc = sqr(lc)*(sqr(la) + sqr(lb) - sqr(lc));
  real sum = bca + bcb + bcc;
  if (fabs(sum) < 1e-6) {
    x = (1.0/3.0)*(a + b + c);
    radius = 1e99;
  } else {
    x = bca*a + bcb*b + bcc*c;
    x *= 1.0/sum;
    radius = (x-a).abs();
  }
}

vec3_t getBarycentricCoordinates(real x, real y)
{
  if(isnan(x) || isinf(x) || isnan(y) || isinf(y)) {
    cout << "x="<<x;
    cout << "y="<<y;
    BUG;
  }
  
  real x_1=0;
  real y_1=0;
  real x_2=1;
  real y_2=0;
  real x_3=0;
  real y_3=1;
  
  mat2_t T;
  T[0][0]=x_1-x_3; T[0][1]=x_2-x_3;
  T[1][0]=y_1-y_3; T[1][1]=y_2-y_3;
  
  if(T.det()==0) {
    cout << "T.det()="<<T.det();
    cout << T[0][0]<<T[0][1];
    cout << T[1][0]<<T[1][1];
    cout << "T[0][0]*T[1][1]-T[1][0]*T[0][1]="<<T[0][0]*T[1][1]-T[1][0]*T[0][1];
    BUG;
  }
  
  real lambda_1 = ((y_2-y_3)*(x-x_3)-(x_2-x_3)*(y-y_3))/(T.det());
  real lambda_2 = (-(y_1-y_3)*(x-x_3)+(x_1-x_3)*(y-y_3))/(T.det());
  real lambda_3 = 1-lambda_1-lambda_2;
  
  vec3_t bary_coords(lambda_1,lambda_2,lambda_3);
  return bary_coords;
  
}

vec3_t intersectionOnPlane(vec3_t v, vec3_t A, vec3_t nA, vec3_t B, vec3_t nB)
{
  vec3_t u = B-A;
  v.normalise();
  v = u.abs()*v;
  
  vec2_t p_A(0,0);
  vec2_t p_B(1,0);
  vec2_t p_nA = projectVectorOnPlane(nA,u,v);
  vec2_t p_nB = projectVectorOnPlane(nB,u,v);
  
  vec2_t p_tA = turnRight(p_nA);
  vec2_t p_tB = turnRight(p_nB);
  
  real k1, k2;
  vec2_t p_K;
  if (!intersection(k1, k2, p_A, p_tA, p_B, p_tB)) {
    p_K = 0.5*(p_A + p_B);
  } else {
    p_K = p_A + k1*p_tA;
  }
  
  if (p_K[0] <0 ) {
    p_K[0] = 0;
  }
  if (p_K[0] >1 ) {
    p_K[0] = 1;
  }
  vec3_t K = A + p_K[0]*u + p_K[1]*v;
  return K;
}

vec2_t projectVectorOnPlane(vec3_t V,vec3_t i,vec3_t j)
{
  if (i.abs2()==0) BUG;
  if (j.abs2()==0) BUG;
  real x = V*i/i.abs2();
  real y = V*j/j.abs2();
  return vec2_t(x,y);
}

vec3_t projectPointOnPlane(const vec3_t& M, const vec3_t& A, const vec3_t& N)
{
  real k = ((M-A)*N)/N.abs2();
  return( M - k*N );
}

vec3_t projectPointOnEdge(const vec3_t& M,const vec3_t& A, const vec3_t& u)
{
  checkVector(u);
  if (u.abs2()==0) BUG;
  real k = ((M-A)*u)/u.abs2();
  return A + k*u;
}

void cart2spherical(vec3_t x, real &alpha, real &beta, real &r)
{
  r = x.abs();
  static const vec3_t ex(1,0,0);
  vec3_t xy(x[0],x[1],0);
  if (x[1] >= 0) {
    alpha = angle(ex, xy);
  } else {
    alpha = 2*M_PI - angle(xy, ex);
  }
  if (xy.abs2() > 0) {
    if (x[2] >= 0) {
      beta = angle(xy, x);
    } else {
      beta = -angle(xy, x);
    }
  } else {
    beta = 0.5*M_PI;
  }
}

bool isInsideCartesianBox(vec3_t x, vec3_t x1, vec3_t x2)
{
  for (int i = 0; i < 3; ++i) {
    if (x[i] < x1[i] || x[i] > x2[i]) {
      return false;
    }
  }
  return true;
}

} // namespace
