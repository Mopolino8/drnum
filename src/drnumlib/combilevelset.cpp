#include "combilevelset.h"


CombiLevelSet::CombiLevelSet(LevelSetDefinition* levelset_a, LevelSetDefinition* levelset_b)
{
  m_LevelSets.clear();
  m_LevelSets.push_back(levelset_a);
  m_LevelSets.push_back(levelset_b);
  findLowestLevelSets();

  //  m_LevelSetA = levelset_a;
  //  m_LevelSetB = levelset_b;
  //  findLowestLevelSets();
}


CombiLevelSet::CombiLevelSet(LevelSetDefinition* levelset_a)
{
  m_LevelSets.clear();
  m_LevelSets.push_back(levelset_a);
  findLowestLevelSets();
}


void CombiLevelSet::includeLevelSet(LevelSetDefinition* levelset)
{
  m_LevelSets.push_back(levelset);
  considerLowestLevelSetsOf(levelset);
}


void CombiLevelSet::considerLowestLevelSetsOf(LevelSetDefinition* levelset)
{
  vector<LevelSetDefinition*> your_lowest_levelsets;
  levelset->getLowestLevelSets(your_lowest_levelsets);
  concatLowestLevelSets(your_lowest_levelsets);
}


void CombiLevelSet::concatLowestLevelSets(vector<LevelSetDefinition*>& other_lowest_levelsets)
{
  // Append elements of other_lowest_levelsets onto own list, then sort and make unique
  //.. append
  for (size_t i2 = 0; i2 < other_lowest_levelsets.size(); i2++) {
    m_LowestLevelSets.push_back(other_lowest_levelsets[i2]);
  }
  //.. sort
  sort(m_LowestLevelSets.begin(), m_LowestLevelSets.end());
  //.. remove duplicates and resize
  vector<LevelSetDefinition*>::iterator it;
  it = unique(m_LowestLevelSets.begin(), m_LowestLevelSets.end());
  m_LowestLevelSets.resize(it - m_LowestLevelSets.begin());
}


void CombiLevelSet::findLowestLevelSets()
{
  m_LowestLevelSets.clear();

  // Loop for LevelSets in this combo.
  // Let levelsets write into help list and concatenate on m_LowestLevelSets
  for (size_t i_o = 0; i_o < m_LevelSets.size(); i_o++) {
    vector<LevelSetDefinition*> your_lowest_levelsets;
    m_LevelSets[i_o]->getLowestLevelSets(your_lowest_levelsets);
    concatLowestLevelSets(your_lowest_levelsets);
  }

//  //.. Let levelset A write directly into own list of "this"
//  m_LevelSetA->getLowestLevelSets(m_LowestLevelSets);

//  //.. Let levelset B write into help list and concatenate on m_LowestLevelSets
//  vector<LevelSetDefinition*> your_lowest_levelsets;
//  m_LevelSetB->getLowestLevelSets(your_lowest_levelsets);
//  concatLowestLevelSets(your_lowest_levelsets);
}


real CombiLevelSet::calcDistance(const real& xo, const real& yo, const real& zo)
{
  // Check wether point is inside lowest levelsets
  for (size_t i_low = 0; i_low < m_LowestLevelSets.size(); i_low++) {
    real low_dist = m_LowestLevelSets[i_low]->calcDistance(xo, yo, zo);
    m_LowestLevelSets[i_low]->setKnownDistance(low_dist);
  }

  // Own operation
  real distance = evalReal();

  return distance;
}
