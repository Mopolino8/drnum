#ifndef COMBILEVELSETAND_H
#define COMBILEVELSETAND_H

class CombiLevelSetAnd;

#include "combilevelset.h"

/**
  * Combine two or more objects in "AND" definition.
  *
  * Geometry:
  * Combined object occupies space commonly covered by all of the input objects.
  *
  * The combined levelset value (the distance from combined object) is the highest
  * of all levelset objects to consider:
  *
  *      G = MAX(G_a, B_bi) ;  all i from the set of lsobjs to consider
  */
class CombiLevelSetAnd : public CombiLevelSet
{
public:
  CombiLevelSetAnd(LevelSetDefinition* levelset_a, LevelSetDefinition* levelset_b);
  CombiLevelSetAnd(LevelSetDefinition* levelset_a);
  virtual real evalReal();
};

#endif // COMBILEVELSETAND_H
