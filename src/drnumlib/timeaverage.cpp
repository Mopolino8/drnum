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

#include "timeaverage.h"

TimeAverage::TimeAverage(real dt, real t_start)
{
  m_Valid = false;
  m_TimeInterval = dt;
  m_AverageValue = 0;
  m_StartTime = t_start;
  m_LastTime = t_start;
  m_TotalValue = 0;
}

void TimeAverage::add(real t, real v)
{
  if (m_TimeInterval > 0) {
    entry_t new_entry;
    new_entry.t = t;
    new_entry.v = v;
    m_Values.append(new_entry);
    if (m_Values.size() > 2) {
      if (m_Values.last().t - m_Values[1].t >= m_TimeInterval) {
        m_Values.removeFirst();
      }
      real T = m_Values.last().t - m_Values.first().t;
      if (T >= m_TimeInterval) {
        real dt = m_TimeInterval - T;
        real w = dt/(m_Values[1].t - m_Values[0].t);
        real t0 = w*m_Values[1].t + (1.0-w)*m_Values[0].t;
        real v0 = w*m_Values[1].v + (1.0-w)*m_Values[0].v;
        QList<entry_t>::iterator i = m_Values.begin();
        ++i;
        real v_integral = 0;
        while (i != m_Values.end()) {
          real t = i->t;
          real v = i->v;
          v_integral += 0.5*(v + v0)*(t - t0);
          t0 = t;
          v0 = v;
          ++i;
        }
        m_AverageValue = v_integral/m_TimeInterval;
        m_Valid = true;
      }
    }
  } else {
    m_TotalValue += (t - m_LastTime)*v;
    if (t - m_StartTime > 0) {
      m_AverageValue = m_TotalValue/(t - m_StartTime);
      m_Valid = true;
    }
  }
  m_LastTime = t;
}
