#ifndef COMPRESSIBLEOBJECT_H
#define COMPRESSIBLEOBJECT_H

class CompressibleObject
{

protected: // attributes

  /*
  restrict real *m_Rho;
  restrict real *m_Rhou;
  restrict real *m_Rhov;
  restrict real *m_Rhow;
  restrict real *m_RhoE;
  restrict real *m_ResRho;
  restrict real *m_ResRhou;
  restrict real *m_ResRhov;
  restrict real *m_ResRhow;
  restrict real *m_ResRhoE;
  */

public:

  real gasR()     { return 287; }    ///< @todo find a concept for this
  real gasGamma() { return 1.4; }    ///< @todo find a concept for this
  real gasCp()    { return 1004.5; } ///< @todo find a concept for this
  real gasCv()    { return 717.5; }  ///< @todo find a concept for this

};

#endif // COMPRESSIBLEOBJECT_H
