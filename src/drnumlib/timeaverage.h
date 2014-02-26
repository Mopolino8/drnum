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
// + DrNUM is distributed in the hope that it will be useful,             +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with DrNUM. If not, see <http://www.gnu.org/licenses/>.        +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef TIMEAVERAGE_H
#define TIMEAVERAGE_H

#include "drnum.h"
#include <QList>

class TimeAverage
{

private: // attributes

  struct entry_t
  {
    real t, v;
  };

  QList<entry_t> m_Values;
  real           m_AverageValue;
  bool           m_Valid;
  real           m_TimeInterval;
  real           m_StartTime;
  real           m_LastTime;
  real           m_TotalValue;

public:

  TimeAverage(real dt = -1.0, real t_start = 0.0);

  void add(real t, real v);
  void setTimeInterval(real dt) { m_TimeInterval = dt; }
  real time() { return m_LastTime; }
  real value() { return m_AverageValue; }
  bool valid() { return m_Valid; }

};

#endif // TIMEAVERAGE_H
