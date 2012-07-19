#ifndef COMPRESSIBLEOBJECT_H
#define COMPRESSIBLEOBJECT_H

class CompressibleObject
{

public:

  struct flux_t
  {
    real rho;
    real rhou;
    real rhov;
    real rhow;
    real rhoE;
  };

  real gasR()     { return 287; }    ///< @todo find a concept for this
  real gasGamma() { return 1.4; }    ///< @todo find a concept for this
  real gasCp()    { return 1004.5; } ///< @todo find a concept for this
  real gasCv()    { return 717.5; }  ///< @todo find a concept for this

};

#endif // COMPRESSIBLEOBJECT_H
