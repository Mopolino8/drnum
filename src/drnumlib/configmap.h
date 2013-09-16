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
#ifndef CONFIGMAP_H
#define CONFIGMAP_H

#include <map>
#include <string>

#include "blockcfd.h"
#include "stringtools.h"


class ConfigMap
{

  std::map<string, string> m_RawMap;

public:

  ConfigMap();

  void readFromFile(string file_name);

  template <typename T> void getValue(string key, T& value);

};


template <typename T>
inline void ConfigMap::getValue(string key, T& value)
{
  StringTools::stringTo(key, value);
}

#endif // CONFIGMAP_H
