#ifndef ROE_H
#define ROE_H

#include "fluxes/compressibleflux.h"
#include "cartesianpatch.h"

template <typename TReconstruction, typename TGas>
class Roe : public CompressibleFlux
{

  TReconstruction* m_Reconstruction;

public: // methods

  Roe(TReconstruction* reconstruction) { m_Reconstruction = reconstruction; }

  real correctLambda(real L);

  void xField(CartesianPatch *patch,
              size_t i, size_t j, size_t k,
              real x, real y, real z,
              real A, real* flux);
  void yField(CartesianPatch *patch,
              size_t i, size_t j, size_t k,
              real x, real y, real z,
              real A, real* flux);
  void zField(CartesianPatch *patch,
              size_t i, size_t j, size_t k,
              real x, real y, real z,
              real A, real* flux);


};


template <typename TReconstruction, typename TGas>
real Roe<TReconstruction, TGas>::correctLambda(real L)
{
  static const real delta = 0.1;
  L = fabs(L);
  if (L < delta) {
    countFlops(3);
    return 0.5*(sqr(L)/delta + delta);
  }
  return L;
}

template <typename TReconstruction, typename TGas>
void Roe<TReconstruction, TGas>::xField(CartesianPatch *patch,
                                        size_t i, size_t j, size_t k,
                                        real x, real y, real z,
                                        real A, real* flux)
{
  COMPRESSIBLE_LEFT_PROJX;
  COMPRESSIBLE_RIGHT_PROJX;

  real r = sqrt(r_l*r_r);
  real D = sqrt(r_r/r_l);
  real u = (D*u_r + u_l)/(D + 1);
  //real u = 0.5*(u_r + u_l);
  real v = 0.5*(v_r + v_l);
  real w = 0.5*(w_r + w_l);
  real H = (D*H_r + H_l)/(D + 1);
  //real H = 0.5*(H_r + H_l);
  real h  = H - 0.5*(sqr(u) + sqr(v) + sqr(w));
  real p = 0.5*(p_r + p_l);
  countSqrts(2);
  countFlops(20);

  real var[5] = {r, r*u, r*v, r*w, r*H - p};
  real a  = sqrt((TGas::gamma(var)-1)*h);
  real a1 = -1/sqr(a);
  real a2 = 0.5/sqr(a);
  countSqrts(1);
  countFlops(10);

  real R0[5] = { a1         , a2           , 0   , 0   , a2           };
  real R1[5] = { a1*u       , a2*(u + a)   , 0   , 0   , a2*(u - a)   };
  real R2[5] = { a1*v       , a2*v         , r   , r   , a2*v         };
  real R3[5] = { a1*w       , a2*w         , r   , r   , a2*w         };
  real R4[5] = { 0.5*a1*u*u , a2*(H + a*u) , r*v , r*w , a2*(H - a*u) };
  countFlops(22);

  real Dr = r_r - r_l;
  real Du = u_r - u_l;
  real Dv = v_r - v_l;
  real Dw = w_r - w_l;
  real Dp = p_r - p_l;
  countFlops(5);

  real DW0 = Dp - a*a*Dr;
  real DW1 = Dp + r*a*Du;
  real DW2 = Dv;
  real DW3 = Dw;
  real DW4 = Dp - r*a*Du;
  countFlops(9);

  real F0_l = ru_l;
  real F1_l = u_l*ru_l + p_l;
  real F2_l = u_l*rv_l;
  real F3_l = u_l*rw_l;
  real F4_l = u_l*(rE_l + p_l);
  countFlops(6);

  real F0_r = ru_r;
  real F1_r = u_r*ru_r + p_r;
  real F2_r = u_r*rv_r;
  real F3_r = u_r*rw_r;
  real F4_r = u_r*(rE_r + p_r);
  countFlops(6);

  real L0 = correctLambda(u);
  real L1 = correctLambda(u + a);
  real L2 = L0;
  real L3 = L0;
  real L4 = correctLambda(u - a);
  countFlops(6);

  flux[0] += A*(0.5*(F0_l + F0_r) - 0.5*(R0[0]*L0*DW0 + R0[1]*L1*DW1 + R0[2]*L2*DW2 + R0[3]*L3*DW3 + R0[4]*L4*DW4));
  flux[1] += A*(0.5*(F1_l + F1_r) - 0.5*(R1[0]*L0*DW0 + R1[1]*L1*DW1 + R1[2]*L2*DW2 + R1[3]*L3*DW3 + R1[4]*L4*DW4));
  flux[2] += A*(0.5*(F2_l + F2_r) - 0.5*(R2[0]*L0*DW0 + R2[1]*L1*DW1 + R2[2]*L2*DW2 + R2[3]*L3*DW3 + R2[4]*L4*DW4));
  flux[3] += A*(0.5*(F3_l + F3_r) - 0.5*(R3[0]*L0*DW0 + R3[1]*L1*DW1 + R3[2]*L2*DW2 + R3[3]*L3*DW3 + R3[4]*L4*DW4));
  flux[4] += A*(0.5*(F4_l + F4_r) - 0.5*(R4[0]*L0*DW0 + R4[1]*L1*DW1 + R4[2]*L2*DW2 + R4[3]*L3*DW3 + R4[4]*L4*DW4));
  countFlops(100);
}

template <typename TReconstruction, typename TGas>
void Roe<TReconstruction, TGas>::yField(CartesianPatch *patch,
                                        size_t i, size_t j, size_t k,
                                        real x, real y, real z,
                                        real A, real* flux)
{
  COMPRESSIBLE_LEFT_PROJY;
  COMPRESSIBLE_RIGHT_PROJY;

  real r = sqrt(r_l*r_r);
  real D = sqrt(r_r/r_l);
  real u = 0.5*(u_r + u_l);
  real v = (D*v_r + v_l)/(D + 1);
  //real v = 0.5*(v_r + v_l);
  real w = 0.5*(w_r + w_l);
  real H = (D*H_r + H_l)/(D + 1);
  //real H = 0.5*(H_r + H_l);
  real h  = H - 0.5*(sqr(u) + sqr(v) + sqr(w));
  real p = 0.5*(p_r + p_l);
  countSqrts(2);
  countFlops(20);

  real var[5] = {r, r*u, r*v, r*w, r*H - p};
  real a  = sqrt((TGas::gamma(var)-1)*h);
  real a1 = -1/sqr(a);
  real a2 = 0.5/sqr(a);
  countSqrts(1);
  countFlops(10);

  real R0[5] = { a1         , 0   , a2           , 0   , a2           };
  real R1[5] = { a1*u       , r   , a2*u         , r   , a2*u         };
  real R2[5] = { a1*v       , 0   , a2*(v + a)   , 0   , a2*(v - a)   };
  real R3[5] = { a1*w       , r   , a2*w         , r   , a2*w         };
  real R4[5] = { 0.5*a1*v*v , r*u , a2*(H + a*v) , r*w , a2*(H - a*v) };
  countFlops(22);

  real Dr = r_r - r_l;
  real Du = u_r - u_l;
  real Dv = v_r - v_l;
  real Dw = w_r - w_l;
  real Dp = p_r - p_l;
  countFlops(5);

  real DW0 = Dp - a*a*Dr;
  real DW1 = Du;
  real DW2 = Dp + r*a*Dv;
  real DW3 = Dw;
  real DW4 = Dp - r*a*Dv;
  countFlops(9);

  real F0_l = rv_l;
  real F1_l = v_l*ru_l;
  real F2_l = v_l*rv_l + p_l;
  real F3_l = v_l*rw_l;
  real F4_l = v_l*(rE_l + p_l);
  countFlops(6);

  real F0_r = rv_r;
  real F1_r = v_r*ru_r;
  real F2_r = v_r*rv_r + p_r;
  real F3_r = v_r*rw_r;
  real F4_r = v_r*(rE_r + p_r);
  countFlops(6);

  real L0 = correctLambda(v);
  real L1 = L0;
  real L2 = correctLambda(v + a);
  real L3 = L0;
  real L4 = correctLambda(v - a);
  countFlops(6);

  flux[0] += A*(0.5*(F0_l + F0_r) - 0.5*(R0[0]*L0*DW0 + R0[1]*L1*DW1 + R0[2]*L2*DW2 + R0[3]*L3*DW3 + R0[4]*L4*DW4));
  flux[1] += A*(0.5*(F1_l + F1_r) - 0.5*(R1[0]*L0*DW0 + R1[1]*L1*DW1 + R1[2]*L2*DW2 + R1[3]*L3*DW3 + R1[4]*L4*DW4));
  flux[2] += A*(0.5*(F2_l + F2_r) - 0.5*(R2[0]*L0*DW0 + R2[1]*L1*DW1 + R2[2]*L2*DW2 + R2[3]*L3*DW3 + R2[4]*L4*DW4));
  flux[3] += A*(0.5*(F3_l + F3_r) - 0.5*(R3[0]*L0*DW0 + R3[1]*L1*DW1 + R3[2]*L2*DW2 + R3[3]*L3*DW3 + R3[4]*L4*DW4));
  flux[4] += A*(0.5*(F4_l + F4_r) - 0.5*(R4[0]*L0*DW0 + R4[1]*L1*DW1 + R4[2]*L2*DW2 + R4[3]*L3*DW3 + R4[4]*L4*DW4));
  countFlops(100);
}

template <typename TReconstruction, typename TGas>
void Roe<TReconstruction, TGas>::zField(CartesianPatch *patch,
                                        size_t i, size_t j, size_t k,
                                        real x, real y, real z,
                                        real A, real* flux)
{
  COMPRESSIBLE_LEFT_PROJZ;
  COMPRESSIBLE_RIGHT_PROJZ;

  real r = sqrt(r_l*r_r);
  real D = sqrt(r_r/r_l);
  real u = 0.5*(u_r + u_l);
  real v = 0.5*(v_r + v_l);
  real w = (D*w_r + w_l)/(D + 1);
  //real w = 0.5*(w_r + w_l);
  real H = (D*H_r + H_l)/(D + 1);
  //real H = (H_r + H_l);
  real h  = H - 0.5*(sqr(u) + sqr(v) + sqr(w));
  real p = 0.5*(p_r + p_l);
  countSqrts(2);
  countFlops(20);

  real var[5] = {r, r*u, r*v, r*w, r*H - p};
  real a  = sqrt((TGas::gamma(var)-1)*h);
  real a1 = -1/sqr(a);
  real a2 = 0.5/sqr(a);
  countSqrts(1);
  countFlops(10);

  real R0[5] = { a1         , 0   , 0   , a2           , a2           };
  real R1[5] = { a1*u       , r   , r   , a2*u         , a2*u         };
  real R2[5] = { a1*v       , r   , r   , a2*v         , a2*v         };
  real R3[5] = { a1*w       , 0   , 0   , a2*(w + a)   , a2*(w - a)   };
  real R4[5] = { 0.5*a1*w*w , r*u , r*v , a2*(H + a*w) , a2*(H - a*w) };
  countFlops(22);

  real Dr = r_r - r_l;
  real Du = u_r - u_l;
  real Dv = v_r - v_l;
  real Dw = w_r - w_l;
  real Dp = p_r - p_l;
  countFlops(5);

  real DW0 = Dp - a*a*Dr;
  real DW1 = Du;
  real DW2 = Dv;
  real DW3 = Dp + r*a*Dw;
  real DW4 = Dp - r*a*Dw;
  countFlops(9);

  real F0_l = rw_l;
  real F1_l = w_l*ru_l;
  real F2_l = w_l*rv_l;
  real F3_l = w_l*rw_l + p_l;
  real F4_l = w_l*(rE_l + p_l);
  countFlops(6);

  real F0_r = rw_r;
  real F1_r = w_r*ru_r;
  real F2_r = w_r*rv_r;
  real F3_r = w_r*rw_r + p_r;
  real F4_r = w_r*(rE_r + p_r);
  countFlops(6);

  real L0 = correctLambda(w);
  real L1 = L0;
  real L2 = L0;
  real L3 = correctLambda(w + a);
  real L4 = correctLambda(w - a);
  countFlops(6);

  flux[0] += A*(0.5*(F0_l + F0_r) - 0.5*(R0[0]*L0*DW0 + R0[1]*L1*DW1 + R0[2]*L2*DW2 + R0[3]*L3*DW3 + R0[4]*L4*DW4));
  flux[1] += A*(0.5*(F1_l + F1_r) - 0.5*(R1[0]*L0*DW0 + R1[1]*L1*DW1 + R1[2]*L2*DW2 + R1[3]*L3*DW3 + R1[4]*L4*DW4));
  flux[2] += A*(0.5*(F2_l + F2_r) - 0.5*(R2[0]*L0*DW0 + R2[1]*L1*DW1 + R2[2]*L2*DW2 + R2[3]*L3*DW3 + R2[4]*L4*DW4));
  flux[3] += A*(0.5*(F3_l + F3_r) - 0.5*(R3[0]*L0*DW0 + R3[1]*L1*DW1 + R3[2]*L2*DW2 + R3[3]*L3*DW3 + R3[4]*L4*DW4));
  flux[4] += A*(0.5*(F4_l + F4_r) - 0.5*(R4[0]*L0*DW0 + R4[1]*L1*DW1 + R4[2]*L2*DW2 + R4[3]*L3*DW3 + R4[4]*L4*DW4));
  countFlops(100);
}

#endif // ROE_H
