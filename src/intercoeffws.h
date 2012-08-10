#ifndef INTERCOEFFWS_H
#define INTERCOEFFWS_H

#include <cstddef>
#include "blockcfd.h"

struct SingleInterCoeffWS;
struct InterCoeffWS;
//class InterCoeffWS;

#include "patch.h"

/// @todo Structures need some better naming conventions to allow easyer understanding

struct SingleInterCoeffWS
{
  size_t indirect_receiveindex;         ///< Index in receive_cells of "this"-owning Patch
  size_t direct_receiveindex;           ///< Cell index in "this"-owning Patch, receiving data
  WeightedSet<real> donor_contribution; ///< WeightedSet containing contributions of donor patch
};

/**
 * Data structure to hold dependency coefficients for data exchange from ONE neighbouring donor patch.
 *
 * Parallel/Vector info: - Unique size_t-sequence: no write conflicts using
 *                         m_DonorCellContribsWS[...].direct_receiveindex
 *                         for left side addressing. Force vectorization, as needed.
 *                       - Potential recurrence for different donor patches (several instances
 *                         of InterCoeffWS for same "this"-owning (receiving) patch). Same cells may be
 *                         receiving data from more than one giving Patch.
 */
struct InterCoeffWS
//class InterCoeffWS
{
//public:
  /**
   * Patch from which data will be received.
   */
  Patch* m_DonorPatch;

//  /**
//   * Cells in "this"-owning Patch, receiving data from giving Patch.
//   * Parallel/Vector info: * Unique size_t-arrays: no write conflicts using m_ReceivingCells
//   *                         for left side array addressing. Force vectorization, as needed.
//   *                       * Potential recurrence for different giving Patches (several instances
//   *                         of InterCoeff for same receiving Patch). Same cells may be receiving data
//   *                         from more than one giving Patch.
//   */
//  TList<size_t>* m_ReceivingCells;

  /**
   * Contributing pairs(cells, weights) in giving Patch.
   */
  vector<SingleInterCoeffWS> m_DonorCellContribsWS;

  /**
    * Set donor patch.
    * @param donor_patch
    */
  void setDonorPatch(Patch* donor_patch)
  {
    m_DonorPatch = donor_patch;  /// @todo this is dangerous, since it might be overwritten
  }

  /**
    * Insert a WeightedSet for a given indexed request
    * @param i_indirect the indirect index of the requesting cell (index in list of requesting cells for owning patch)
    * @param i_direct the direct index of the requesting cell
    * @param contribution WeightedSet with coefficients related to m_DonorPatch
    */
  void push(const size_t& i_indirect, const size_t& i_direct, const WeightedSet<real>& contribution)
  {
    SingleInterCoeffWS sic_h;
    sic_h.indirect_receiveindex = i_indirect;
    sic_h.direct_receiveindex = i_direct;
    sic_h.donor_contribution = contribution;
    m_DonorCellContribsWS.push_back(sic_h);  /// @todo strange debug error on ddd: cannot access operator[] unless the folowing line is active
    //cout << m_DonorCellContribsWS[0].indirect_receiveindex << endl;
    cout << m_DonorCellContribsWS[m_DonorCellContribsWS.size()-1].indirect_receiveindex;
    cout << " w_set-size = " << m_DonorCellContribsWS[m_DonorCellContribsWS.size()-1].donor_contribution.v.size() << endl;
  }

  /**
    * Prepare for averiging: divide contributions for a receiving cell by number of contributing patches
    * NOTE: A receiving cell might be located in the core region of several donor patches.
    * @param hits a vector containing the number of contributing donor patches.
    */
  void adjust2Average(vector<size_t> hits)
  {
    for(size_t ll=0; ll<m_DonorCellContribsWS.size(); ll++) {
      // indirect address (equal to index in hits)
      size_t ind_rec_h = m_DonorCellContribsWS[ll].indirect_receiveindex;
      if(hits[ind_rec_h] > 1) {
        m_DonorCellContribsWS[ll].donor_contribution *= (1./hits[ind_rec_h]);
      }
    }
  }


};
#endif // INTERCOEFFWS_H

