// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of DrNUM.                                          +
// +                                                                      +
// + Copyright 2013 numrax GmbH, enGits GmbH                              +
// +                                                                      +
// + DrNUM is free software: you can redistribute it and/or modify        +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
  real   axx;                        ///< xx component of transformation matrix from donor to receiver
  real   axy;                        ///< xy component of transformation matrix from donor to receiver
  real   axz;                        ///< xz component of transformation matrix from donor to receiver
  real   ayx;                        ///< yx component of transformation matrix from donor to receiver
  real   ayy;                        ///< yy component of transformation matrix from donor to receiver
  real   ayz;                        ///< yz component of transformation matrix from donor to receiver
  real   azx;                        ///< zx component of transformation matrix from donor to receiver
  real   azy;                        ///< zy component of transformation matrix from donor to receiver
  real   azz;                        ///< zz component of transformation matrix from donor to receiver
};

#endif // DONOR_T_H
