#ifndef INTERCOEFFWS_H
#define INTERCOEFFWS_H

#include <cstddef>
#include "blockcfd.h"
#include "math/coordtransform.h"

struct SingleInterCoeffWS;
class InterCoeffWS;

#include "patch.h"
#include "utility/usparseweightedset.h"

/// @todo Structures need some better naming conventions to allow easyer understanding

struct SingleInterCoeffWS
{
  size_t indirect_receiveindex;         ///< Index in indirect array receive_cells of receiving Patch
  size_t direct_receiveindex;           ///< Cell index in receiving Patch
  USparseWeightedSet<real> donor_contribution; ///< WeightedSet containing contributions of donor patch
  /** @todo Storing indirect_receiveindex might be ommitted, if hit counters (Patch::m_receive_cell_data_hits)
    * work on true addresses. This wastes some memory on Patch classes, but might be easyer to understand.
    */
};

/**
 * Data structure to hold dependency coefficients for data exchange from ONE neighbouring donor patch.
 *
 * Receiving patch: keeps one instance of this class per donor neighbour.
 *
 * Parallel/Vector info: - Unique size_t-sequence: no write conflicts using
 *                         m_DonorCellContribsWS[...].direct_receiveindex
 *                         for left side addressing. Force vectorization, as needed.
 *                       - Potential recurrence for different donor patches (several instances
 *                         of InterCoeffWS for same "this"-owning (receiving) patch). Same cells may be
 *                         receiving data from more than one giving Patch.
 */
//struct InterCoeffWS
class InterCoeffWS
{

  friend class InterCoeffPad;

protected: // attributes

  /// @todo new mem-structure: m_DonorPatch probably obsolete due to borrowed pointers anyway
  /// Patch from which data will be received.
  Patch* m_DonorPatch;

  /// Coordinate transformation from m_DonorPatch to receiving Patch.
  CoordTransform m_ct;

  /// Contributing pairs(cells, weights) in giving Patch.
  vector<SingleInterCoeffWS> m_DonorCellContribsWS;


public: // methods

  /**
    * Set donor patch.
    * @param donor_patch
    */
  void setDonorPatch(Patch* donor_patch);

  /**
    * Set relative coordinate transformation. Required for later access when transfering vectorial data.
    */
  void setCoordTransform(const CoordTransform& ct_donor2owner);

  /**
    * Get number of receiving cells.
    */
  size_t getNumRecCells() const;

  /**
    * Get maximum number of donor cells contributing to any single receiver cell.
    */
  size_t getMaxContrib() const;

  /**
    * Get direct receiving cell index in receiving patch
    * @param ll_rec index in list of contributions
    */
  size_t getDirectRecCell(const size_t& ll_rec) const;

  /**
    * Get WeightedSet-ptr for contribution to one receiving cell
    * @param ll_rec index in list of contributions
    */
  USparseWeightedSet<real>* getWS(const size_t& ll_rec);

  /**
    * Insert a WeightedSet for a given indexed request
    * @param i_indirect the indirect index of the requesting cell (index in list of requesting cells for owning patch)
    * @param i_direct the direct index of the requesting cell
    * @param contribution WeightedSet with coefficients related to m_DonorPatch
    */
  void push(const size_t& i_indirect, const size_t& i_direct, const WeightedSet<real>& contribution);

  /**
    * Prepare for averiging: divide contributions for a receiving cell by number of contributing patches
    * NOTE: A receiving cell might be located in the core region of several donor patches.
    * @param hits a vector containing the number of contributing donor patches.
    */
  void adjust2Average(const vector<size_t>& hits);

  /**
    * Apply WeightedSet pattern to donor_data and add to receive_data.
    * No turning of any types of variables, suitable for scalar variables only.
    * NOTE: donor_data and receive_data are small vectors with one entry per variable type to transfer
    * @param donor_data vector of real pointers to variable data of donor patch
    * @param receiver_data vector of real pointers to variable data of receiving patch
    */
  void transFromTo(const vector<real*>& donor_data, const vector<real*>& receiver_data);

  /**
    * Apply WeightedSet pattern to donor_data and add to receive_data.
    * Includes turning of vectorial variable (example: speed vector) given by indicees (i_vx, i_vy, i_vz) in donor_data
    * NOTE: donor_data and receive_data are small vectors with one entry per variable type to transfer
    *
    * @todo must improve variable selection for transformFree operation to allow better optimization.
    *
    * @param donor_data vector of real pointers to variable data of donor patch
    * @param receiver_data vector of real pointers to variable data of receiving patch
    * @param i_vx variable index forming the x-comp. of a vectorial variable. Example: "1" for speed (p,u,v,w,T)
    * @param i_vy variable index forming the y-comp. of a vectorial variable. Example: "2" for speed (p,u,v,w,T)
    * @param i_vz variable index forming the z-comp. of a vectorial variable. Example: "3" for speed (p,u,v,w,T)

    */
  void transTurnFromTo(const vector<real*>& donor_data, const vector<real*>& receiver_data,
                       const size_t& i_vx, const size_t& i_vy, const size_t& i_vz);

  /**
    * Apply WeightedSet pattern to donor_data and add to receive_data.
    * Includes turning of vectorial variable (example: speed vector) given by fixed indicees (0, 1, 2) in donor_data
    * NOTE: donor_data and receive_data are small vectors with one entry per variable type to transfer
    * @param donor_data vector of real pointers to variable data of donor patch
    * @param receiver_data vector of real pointers to variable data of receiving patch
    */
  void transTurnFromTo(const vector<real*>& donor_data, const vector<real*>& receiver_data);

};  // end class definition


inline void InterCoeffWS::setDonorPatch(Patch* donor_patch)
{
  m_DonorPatch = donor_patch;
}


inline void InterCoeffWS::setCoordTransform(const CoordTransform& ct_donor2owner)
{
  m_ct = ct_donor2owner;
}

inline size_t InterCoeffWS::getNumRecCells() const
{
  return m_DonorCellContribsWS.size();
}

inline size_t InterCoeffWS::getMaxContrib() const
{
  size_t max_contrib_size = 0;
  for (size_t ll_rec = 0; ll_rec < m_DonorCellContribsWS.size(); ll_rec++) {
    size_t s_h = m_DonorCellContribsWS[ll_rec].donor_contribution.getSize();
    if(max_contrib_size < s_h) {
      max_contrib_size = s_h;
    }
  }
  return max_contrib_size;
}

inline size_t InterCoeffWS::getDirectRecCell(const size_t& ll_rec) const
{
  return m_DonorCellContribsWS[ll_rec].direct_receiveindex;
}

inline USparseWeightedSet<real>* InterCoeffWS::getWS(const size_t& ll_rec)
{
  return &(m_DonorCellContribsWS[ll_rec].donor_contribution);
}

inline void InterCoeffWS::push(const size_t& i_indirect, const size_t& i_direct, const WeightedSet<real>& contribution)
{
  SingleInterCoeffWS sic_h;
  sic_h.indirect_receiveindex = i_indirect;
  sic_h.direct_receiveindex = i_direct;
  sic_h.donor_contribution = contribution;
  m_DonorCellContribsWS.push_back(sic_h);  /// @todo strange debug error on ddd: cannot access operator[] unless the folowing line is active
  // cout << m_DonorCellContribsWS[m_DonorCellContribsWS.size()-1].indirect_receiveindex;
  // cout << " w_set-size = " << m_DonorCellContribsWS[m_DonorCellContribsWS.size()-1].donor_contribution.v.size() << endl;
}


inline void InterCoeffWS::adjust2Average(const vector<size_t>& hits)
{
  for(size_t ll=0; ll<m_DonorCellContribsWS.size(); ll++) {
    // indirect address (equal to index in hits)
    size_t ind_rec_h = m_DonorCellContribsWS[ll].indirect_receiveindex;
    if(hits[ind_rec_h] > 1) {
      m_DonorCellContribsWS[ll].donor_contribution *= (1./hits[ind_rec_h]);
    }
  }
}


inline void InterCoeffWS::transFromTo(const vector<real*>& donor_data, const vector<real*>& receiver_data) {
#ifdef DEBUG
  // check sizes
  if(donor_data.size() != receiver_data.size()) {BUG;}
#endif
  size_t num_vars = donor_data.size();
  for (size_t i_ind = 0; i_ind < m_DonorCellContribsWS.size(); i_ind++) {
    size_t i_receive = m_DonorCellContribsWS[i_ind].direct_receiveindex;
    USparseWeightedSet<real> &dccws = m_DonorCellContribsWS[i_ind].donor_contribution;
    for(size_t i_v=0; i_v<num_vars; i_v++) {
      *(receiver_data[i_v]+i_receive) += dccws.computeValue(donor_data[i_v]);
    }
  }
}


inline void InterCoeffWS::transTurnFromTo(const vector<real*>& donor_data, const vector<real*>& receiver_data,
                                          const size_t& i_vx, const size_t& i_vy, const size_t& i_vz)
{
#ifdef DEBUG
  // check vector indexing bounds
  if(i_vx > donor_data.size()) {BUG;}
  if(i_vy > donor_data.size()) {BUG;}
  if(i_vz > donor_data.size()) {BUG;}
  // check sizes
  if(donor_data.size() != receiver_data.size()) {BUG;}
#endif
  size_t num_vars = donor_data.size();
  real* vars_h = new real[num_vars];
  for(size_t i_ind=0; i_ind<m_DonorCellContribsWS.size(); i_ind++) {
    size_t i_receive = m_DonorCellContribsWS[i_ind].direct_receiveindex;
    USparseWeightedSet<real> &dccws = m_DonorCellContribsWS[i_ind].donor_contribution;
    for(size_t i_v=0; i_v<num_vars; i_v++) {
      vars_h[i_v] = dccws.computeValue(donor_data[i_v]);
    }
    // turn vectorial data onto system of receiver
    /** @todo likely to be a performance issue: fixed indexing for vectorial quantities prefereable, if possible */
    m_ct.transfree(vars_h[i_vx], vars_h[i_vy], vars_h[i_vz]);
    // contribute
    for(size_t i_v=0; i_v<num_vars; i_v++) {
      *(receiver_data[i_v]+i_receive) += vars_h[i_v];
    }
  }
  delete vars_h;
}

inline void InterCoeffWS::transTurnFromTo(const vector<real*>& donor_data, const vector<real*>& receiver_data)
{
#ifdef DEBUG
  // check number of variable types
  if(donor_data.size() != receiver_data.size()) {BUG;}
  if(donor_data.size() < 3) {BUG;}
#endif
  size_t num_vars = donor_data.size();
  for(size_t i_ind=0; i_ind<m_DonorCellContribsWS.size(); i_ind++) {
    size_t i_receive = m_DonorCellContribsWS[i_ind].direct_receiveindex;
    USparseWeightedSet<real> &dccws = m_DonorCellContribsWS[i_ind].donor_contribution;
    // do vectorial var first
    real vec_vars_0 = dccws.computeValue(donor_data[0]);    // x_comp vectorial var to turn
    real vec_vars_1 = dccws.computeValue(donor_data[1]);    // y_comp vectorial var to turn
    real vec_vars_2 = dccws.computeValue(donor_data[2]);    // z_comp vectorial var to turn
    m_ct.transfree(vec_vars_0, vec_vars_1, vec_vars_2);  // turn into system of receiver
    *(receiver_data[0]+i_receive) += vec_vars_0;         // add contribution
    *(receiver_data[1]+i_receive) += vec_vars_1;         //        "
    *(receiver_data[2]+i_receive) += vec_vars_2;         //        "
    // do all other variables
    for(size_t i_v=3; i_v<num_vars; i_v++) {
      *(receiver_data[i_v]+i_receive) += dccws.computeValue(donor_data[i_v]);
    }
  }
}


#endif // INTERCOEFFWS_H

