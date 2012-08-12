#ifndef INTERCOEFFPAD_H
#define INTERCOEFFPAD_H

#include <cstddef>
#include "blockcfd.h"
#include "math/coordtransform.h"
#include "intercoeffws.h"

class InterCoeffPad;

#include "patch.h"


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
class InterCoeffPad
{

protected: // attributes

  /// @todo new mem-structure: m_DonorPatch probably obsolete due to borrowed pointers anyway
  /// Patch from which data will be received.
  Patch* m_DonorPatch;

  /// Coordinate transformation from m_DonorPatch to receiving Patch.
  CoordTransform m_ct;

  /// Number of receiving cells.
  size_t m_NumRecCells;

  /// Numbers of contributing pairs(cells, weights) per receiving cell (constant stride).
  size_t m_StrideGivePerRec;


  /**
   * Cells in "this"-owning Patch, receiving data from donor Patch.
   * Parallel/Vector info: * Unique size_t-arrays: no write conflicts using m_ReceivingCells
   *                         for left side array addressing. Force vectorization, as needed.
   *                       * Potential recurrence for different giving Patches (several instances
   *                         of InterCoeff for same receiving Patch). Same cells may be receiving data
   *                         from more than one giving Patch.
   */
  size_t* m_RecCells;

  /**
   * Contributing cells of donor Patch.
   * Note: Fixed stride. To allow more performant parallel/vector execution, for each receiving cell
   *       in m_ReceivingCells a fixed number m_StrideGivePerRec of pairs(m_DonorCells, m_DonorWeights)
   *       is provided. Padding is used to fill empty spaces, in case of lesser than max contributors.
   */
  size_t* m_DonorCells;

  /**
    * Weights of donor cell contributions. Index match with m_DonorCells.
    * Note: see note for m_DonorCells.
    */
  real* m_DonorWeights;


public: // methods

  /**
    * Constructor.
    * @param icws InterCoeffWS with data content to copy into padded data structure of "this".
    */
  void build(InterCoeffWS& icws);

  /**
    * Apply data pattern to donor_data and add to receive_data.
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


inline void InterCoeffPad::build(InterCoeffWS& icws)
{
  m_DonorPatch = icws.m_DonorPatch;
  m_ct = icws.m_ct;
  // Number of receiving cells
  m_NumRecCells = icws.getNumRecCells();
  // Maximum number of contributors for any receiving cell.
  m_StrideGivePerRec = icws.getMaxContrib();
  // Build up pattern data sets
  m_RecCells = new size_t[m_NumRecCells];
  m_DonorCells = new size_t[m_NumRecCells * m_StrideGivePerRec];
  m_DonorWeights = new real[m_NumRecCells * m_StrideGivePerRec];
  // Fill data into lists m_RecCells, m_DonorCells, m_DonorCells
  size_t start = 0;
  //.. loop through receiving cells
  for (size_t ll_rec = 0; ll_rec < m_NumRecCells; ll_rec++) {
    m_RecCells[ll_rec] = icws.getDirectRecCell(ll_rec); // direct receiving cell index in receiving patch
    WeightedSet<real>* ws_h = icws.getWS(ll_rec);
    for( size_t l_contrib = 0; l_contrib < ws_h->getSize(); l_contrib++) { // fill contents
      m_DonorCells[start+l_contrib] = ws_h->getIndex(l_contrib);
      m_DonorWeights[start+l_contrib] = ws_h->getWeight(l_contrib);
    }
    for( size_t l_pad = ws_h->getSize(); l_pad < m_StrideGivePerRec; l_pad++) { // padding
      m_DonorCells[start+l_pad] = m_DonorCells[start]; // any can do, just be sure not to leave segmentation
      m_DonorWeights[start+l_pad] = 0.;
    }
    start += m_StrideGivePerRec;
  }
}

inline void InterCoeffPad::transFromTo(const vector<real*>& donor_data, const vector<real*>& receiver_data) {
#ifdef DEBUG
  // check sizes
  if(donor_data.size() != receiver_data.size()) {BUG;}
#endif
  // Parallel/Vector info: - No right hand side recurrence (write-conflict) in first loop (ll),
  //                         since m_RecCells is unique.
  //                         Force vector or thread parallel, if required.
  //                       - Identical loop bounds of inner loops (l_contrib, i_v) for all outer loop iterations (ll)
  size_t num_vars = donor_data.size();
  for (size_t ll = 0; ll < m_NumRecCells; ll++) {  // indirect loop for receiving cells
    size_t l_receive = m_RecCells[ll];             // direct receiving cells in receiving patch
    size_t start = ll*m_StrideGivePerRec;          // start address in m_DonorCells/m_DonorWeights pattern
    for(size_t l_contrib = 0; l_contrib < m_StrideGivePerRec; l_contrib++) { // ind loop for contributing cells of donor patch
      size_t donor_cell_index = start + l_contrib;
      size_t donor_cell = m_DonorCells[donor_cell_index];   // direct donor cell index in donor patch
      real donor_weight = m_DonorWeights[donor_cell_index]; // donor cell weight
      for(size_t i_v = 0; i_v < num_vars; i_v++) {          // loop for variables
        *(receiver_data[i_v]+l_receive) += donor_data[i_v][donor_cell] * donor_weight;  // contribute to receiving cell
      }
    }
  }
}


inline void InterCoeffPad::transTurnFromTo(const vector<real*>& donor_data, const vector<real*>& receiver_data,
                                           const size_t& i_vx, const size_t& i_vy, const size_t& i_vz)
{
  BUG;  // not implemented, probably never needed
}


inline void InterCoeffPad::transTurnFromTo(const vector<real*>& donor_data, const vector<real*>& receiver_data)
{
#ifdef DEBUG
  // check number of variable types
  if(donor_data.size() != receiver_data.size()) {BUG;}
  if(donor_data.size() < 3) {BUG;}
#endif
  // Parallel/Vector info: - No right hand side recurrence (write-conflict) in first loop (ll),
  //                         since m_RecCells is unique.
  //                         Force vector or thread parallel, if required.
  //                       - Identical loop bounds of inner loops (l_contrib, i_v) for all outer loop iterations (ll)
  size_t num_vars = donor_data.size();
  for (size_t ll = 0; ll < m_NumRecCells; ll++) {  // indirect loop for receiving cells
    size_t l_receive = m_RecCells[ll];             // direct receiving cells in receiving patch
    size_t start = ll*m_StrideGivePerRec;          // start address in m_DonorCells/m_DonorWeights pattern
    real vec_vars_0 = 0.;    // x_comp vectorial var to turn
    real vec_vars_1 = 0.;    // y_comp vectorial var to turn
    real vec_vars_2 = 0.;    // z_comp vectorial var to turn
    for(size_t l_contrib = 0; l_contrib < m_StrideGivePerRec; l_contrib++) { // ind loop for contributing cells of donor patch
      size_t donor_cell_index = start + l_contrib;
      size_t donor_cell = m_DonorCells[donor_cell_index];   // direct donor cell index in donor patch
      real donor_weight = m_DonorWeights[donor_cell_index]; // donor cell weight
      vec_vars_0 += donor_data[0][donor_cell] * donor_weight;
      vec_vars_1 += donor_data[1][donor_cell] * donor_weight;
      vec_vars_2 += donor_data[2][donor_cell] * donor_weight;
      for(size_t i_v=3; i_v<num_vars; i_v++) { // loop for remaining variables
        *(receiver_data[i_v]+l_receive) += donor_data[i_v][donor_cell] * donor_weight;  // contribute to receiving cell
      }
    }
    m_ct.transfree(vec_vars_0, vec_vars_1, vec_vars_2);  // turn vectorial variables into system of receiver
    *(receiver_data[0]+l_receive) += vec_vars_0;         // add contribution
    *(receiver_data[1]+l_receive) += vec_vars_1;         //        "
    *(receiver_data[2]+l_receive) += vec_vars_2;         //        "
  }
}


#endif // INTERCOEFFPAD_H

