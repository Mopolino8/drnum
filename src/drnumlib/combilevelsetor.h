#ifndef COMBILEVELSETOR_H
#define COMBILEVELSETOR_H

class CombiLevelSetOr;

#include "combilevelset.h"

/**
  * Combine two or more objects in "or" definition.
  *
  * Geometry:
  * Combined object occupies space covered by any of the input objects.
  *
  * The combined levelset value (the distance from combined object) is the lowest
  * of all levelset objects to consider:
  *
  *      G = MIN(G_a, B_bi) ;  all i from the set of lsobjs to consider
  */
class CombiLevelSetOr : public CombiLevelSet
{
public:
    CombiLevelSetOr(LevelSetDefinition* levelset_a, LevelSetDefinition* levelset_b);
    CombiLevelSetOr(LevelSetDefinition* levelset_a);
    virtual real evalReal();
};

#endif // COMBILEVELSETOR_H
