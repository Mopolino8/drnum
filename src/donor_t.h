#ifndef DONOR_T_H
#define DONOR_T_H

/**
  * Struct holding data for one donor->receiver relation.
  * To be owned by receiving patch.
  */
struct donor_t
{
  size_t variable_size;              ///< number of reals per variable in donor patch
  real*  data;                       ///< pointer to m_Data of donor patch (on device being!)
  size_t num_receiver_cells;         ///< Number of cells, receiving data from donor patch
  size_t stride;                     ///< Fixed number of donor cells for each receiving cell
  size_t receiver_index_field_start; ///< Starting index in concatenated receiving cell indicees field of receiving patch
  size_t donor_wi_field_start;       ///< Starting index in concatenated index and weight field for all donor patches
};

#endif // DONOR_T_H
