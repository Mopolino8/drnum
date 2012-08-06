#ifndef INTERCOEFF_H
#define INTERCOEFF_H

#include <cstddef>
#include "blockcfd.h"

//struct InterCoeff;
class InterCoeff;

#include "patch.h"

/**
 * Data structure to hold dependency coefficients for data exchange between patches
 *
 * NOTE: * Dependencies from ONE other Patch only.
 */
struct InterCoeff
    //class InterCoeff
{
  /**
   * Patch from which data will be received.
   */
  Patch* m_GivingPatch;

  /**
   * Cells in "this"-owning Patch, receiving data from giving Patch.
   * Parallel/Vector info: * Unique size_t-arrays: no write conflicts using m_ReceivingCells
   *                         for left side array addressing. Force vectorization, as needed.
   *                       * Potential recurrence for different giving Patches (several instances
   *                         of InterCoeff for same receiving Patch). Same cells may be receiving data
   *                         from more than one giving Patch.
   */
  size_t* m_ReceivingCells;

  /**
   * Numbers of contributing pairs(cells, weights) per receiving cell (constant stride). 
   */
  size_t m_StrideGivePerRec;

  /**
   * Contributing pairs(cells, weights) in giving Patch.
   * Note: Fixed stride. To allow easyer parallel/vector execution, for each receiving cell in
   *       in m_ReceivingCells a fixed number m_StrideGivePerRec of pairs(cells, weights)
   *       is provided. Padding is used to fill empty spaces, in case of lesser contributors.
   */
  pair<size_t, real>* m_GivingCellContribs;

};
#endif // INTERCOEFF_H

