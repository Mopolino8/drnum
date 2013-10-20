#ifndef COMBILEVELSETANDNOT_H
#define COMBILEVELSETANDNOT_H

class CombiLevelSetAndNot;

#include "combilevelset.h"

/**
  * Combine two or more objects in "AND NOT" definition.
  *
  * Geometry:
  * Combined object occupies space covered by object A, unless also covered by
  * any of the objects B_i, all i.
  *
  * The levelset value (the distance) is
  *      G = MAX(G_a, -G_bi) ;  all i from the set to exclude
  */
class CombiLevelSetAndNot : public CombiLevelSet
{
public:
  CombiLevelSetAndNot(LevelSetDefinition* levelset_a, LevelSetDefinition* levelset_b);
  CombiLevelSetAndNot(LevelSetDefinition* levelset_a);
  virtual real evalReal();
};

#endif // COMBILEVELSETANDNOT_H
