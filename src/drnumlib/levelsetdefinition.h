#ifndef LEVELSETDEFINITION_H
#define LEVELSETDEFINITION_H

class LevelSetDefinition;

#include "blockcfd.h"

/**
  * Base class for the definition of levelset objects.
  */
class LevelSetDefinition
{

protected: // attributes
  real m_KnownDistance;

protected: // methods

  /** @todo define wether to use o-coords, or generalized parent coords xyz
    *       in all derived class. */

public:

  void setKnownDistance(const real& dist) {m_KnownDistance = dist;}
  real getKnownDistance() {return m_KnownDistance;}

  LevelSetDefinition();

  virtual real calcDistance(const real& xo, const real& yo, const real& zo) {
    cout << " need derived class for LevelSetDefinition::calcDistance" << endl;
    BUG;
  }

  real calcDistance(vec3_t xyzo) { return calcDistance (xyzo[0], xyzo[1],xyzo[2]); }

  virtual void getLowestLevelSets(vector<LevelSetDefinition*>& my_lowest_levelsets);

  virtual real evalReal() {return getKnownDistance();}

};

#endif // LEVELSETDEFINITION_H
