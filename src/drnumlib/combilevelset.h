#ifndef COMBILEVELSET_H
#define COMBILEVELSET_H

class CombiLevelSet;

#include "levelsetdefinition.h"

class CombiLevelSet : public LevelSetDefinition
{

protected: // attributes
  vector<LevelSetDefinition*> m_LowestLevelSets;
  LevelSetDefinition* m_LevelSetA;
  LevelSetDefinition* m_LevelSetB;

  vector<LevelSetDefinition*> m_LevelSets;

protected: // methods
  void considerLowestLevelSetsOf(LevelSetDefinition* levelset);
  void concatLowestLevelSets(vector<LevelSetDefinition*>& other_lowest_levelsets);
  void findLowestLevelSets();

public:
  CombiLevelSet(LevelSetDefinition* levelset_a, LevelSetDefinition* levelset_b);
  CombiLevelSet(LevelSetDefinition* levelset_a);
  void includeLevelSet(LevelSetDefinition* levelset);
  virtual real calcDistance(const real& xo, const real& yo, const real& zo);
  virtual real evalReal() = 0;

};

#endif // COMBILEVELSET_H
