#ifndef COMPRESSIBLECARTESIANPATCH_H
#define COMPRESSIBLECARTESIANPATCH_H

#include "cartesianpatch.h"
#include "compressibleobject.h"

class CompressibleCartesianPatch : public CartesianPatch, public CompressibleObject
{

protected: // attributes

  //TEulerFlux m_EulerFlux;


protected: // methods

  //void positiveContribution(size_t i, size_t j, size_t k, RealVec<5> flux);
  //void negativeContribution(size_t i, size_t j, size_t k, RealVec<5> flux);


public: // methods

  CompressibleCartesianPatch();

  void setState(size_t i, size_t j, size_t k, real p, real T);
  void getState(size_t i, size_t j, size_t k, real &p, real &T);
  void getState(size_t i, size_t j, size_t k, real &p, real &u, real &v, real& w, real &T);

  //virtual void subStep(real dt);
  //virtual void preStep();


};



/*
template <class TEulerFlux>
inline void CompressibleCartesianPatch<TEulerFlux>::positiveContribution(size_t i, size_t j, size_t k, RealVec<5> flux)
{
  f(i_res, 0, i, j, k) += flux.var[0];
  f(i_res, 1, i, j, k) += flux.var[1];
  f(i_res, 2, i, j, k) += flux.var[2];
  f(i_res, 3, i, j, k) += flux.var[3];
  f(i_res, 4, i, j, k) += flux.var[4];
}

template <class TEulerFlux>
inline void CompressibleCartesianPatch<TEulerFlux>::negativeContribution(size_t i, size_t j, size_t k, RealVec<5> flux)
{
  f(i_res, 0, i, j, k) -= flux.var[0];
  f(i_res, 1, i, j, k) -= flux.var[1];
  f(i_res, 2, i, j, k) -= flux.var[2];
  f(i_res, 3, i, j, k) -= flux.var[3];
  f(i_res, 4, i, j, k) -= flux.var[4];
}
*/

/*
template <class TEulerFlux>
inline void CompressibleCartesianPatch<TEulerFlux>::subStep(real dt)
{
  setFieldToZero(i_res);
  preStep();
  real V = dx()*dy()*dz();
  addField(i_old, dt/V, i_res, 0);
  postStep();
}
*/

/*
template <class TEulerFlux>
inline void CompressibleCartesianPatch<TEulerFlux>::preStep()
{
  m_EulerFlux.setRho (getVariable(0, 0));
  m_EulerFlux.setRhou(getVariable(0, 1));
  m_EulerFlux.setRhov(getVariable(0, 2));
  m_EulerFlux.setRhow(getVariable(0, 3));
  m_EulerFlux.setRhoE(getVariable(0, 4));
  m_EulerFlux.setResRho (getVariable(i_res, 0));
  m_EulerFlux.setResRhou(getVariable(i_res, 1));
  m_EulerFlux.setResRhov(getVariable(i_res, 2));
  m_EulerFlux.setResRhow(getVariable(i_res, 3));
  m_EulerFlux.setResRhoE(getVariable(i_res, 4));

  RealVec<5> flux;

  // x direction
  //
  for (size_t i = 1; i < sizeI(); ++i) {
    for (size_t j = 0; j < sizeJ(); ++j) {
      for (size_t k = 0; k < sizeK(); ++k) {
        if (i == 502) {
          //cout << "break" << endl;
        }
        m_EulerFlux.x(this, i, j, k, dy()*dz(), flux);
        negativeContribution(i-1, j, k, flux);
        positiveContribution(i, j, k, flux);
      }
    }
  }
  for (size_t j = 0; j < sizeJ(); ++j) {
    for (size_t k = 0; k < sizeK(); ++k) {
      real p, T;
      getState(0, j, k, p, T);
      flux.rho  = 0;
      flux.rhou = p*dy()*dz();
      flux.rhov = 0;
      flux.rhow = 0;
      flux.rhoE = 0;
      positiveContribution(0, j, k, flux);
      getState(sizeI() - 1, j, k, p, T);
      flux.rho  = 0;
      flux.rhou = p*dy()*dz();
      flux.rhov = 0;
      flux.rhow = 0;
      flux.rhoE = 0;
      negativeContribution(sizeI() - 1, j, k, flux);
    }
  }

  // y direction
  //
  for (size_t i = 0; i < sizeI(); ++i) {
    for (size_t j = 1; j < sizeJ(); ++j) {
      for (size_t k = 0; k < sizeK(); ++k) {
      }
    }
  }
  for (size_t i = 0; i < sizeI(); ++i) {
    for (size_t k = 0; k < sizeK(); ++k) {
    }
  }

  // z direction
  //
  for (size_t i = 0; i < sizeI(); ++i) {
    for (size_t j = 0; j < sizeJ(); ++j) {
      for (size_t k = 1; k < sizeK(); ++k) {
      }
    }
  }
  for (size_t i = 0; i < sizeI(); ++i) {
    for (size_t j = 0; j < sizeJ(); ++j) {
    }
  }
}
*/

#endif // COMPRESSIBLECARTESIANPATCH_H
