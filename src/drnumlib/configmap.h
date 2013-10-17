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

#include <QMap>
#include <QVariant>

#include "drnum.h"


class ConfigMap
{

private: // attributes

  QMap<QString, QVariant> m_Map;


private: // methods

  void processComments(QString &buffer, QString begin, QString end);
  bool assign(QString key, QString type, QString value);


public:  // methods

  void clear();
  void addFile(QString file_name);
  void addDirectory(QString path);

  template <typename T> void getValue(QString key, T& t)
  {
    QVariant variant = t;
    if (variant.typeName() != m_Map[key].typeName()) {
      QString msg = "Error retrieving parameter \"" + key + "\" of type " + variant.typeName();
      ERROR(qPrintable(msg));
    }
    t = m_Map[key].value<T>();
  }

  template <typename T> T getValue(QString key)
  {
    QVariant variant = T();
    if (variant.typeName() != m_Map[key].typeName()) {
      QString msg = "Error retrieving parameter \"" + key + "\" of type " + variant.typeName();
      ERROR(qPrintable(msg));
    }
    return m_Map[key].value<T>();
  }

  bool exists(QString key);

};


#endif // CONFIGMAP_H
