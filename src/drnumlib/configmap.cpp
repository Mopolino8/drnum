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
#include "configmap.h"

#include <iostream>
#include <fstream>

#include <QDir>
#include <QFile>
#include <QTextStream>
#include <QString>

void ConfigMap::clear()
{
  m_Map.clear();
}

void ConfigMap::processComments(QString &buffer, QString begin, QString end)
{
  QString new_buffer;
  bool comment = false;
  int i = 0;
  while (i < buffer.size()) {
    bool add_to_new_buffer = true;
    if (!comment) {
      if (buffer.mid(i, begin.size()) == begin) {
        comment = true;
        add_to_new_buffer = false;
        i += begin.size() - 1;
      }
    } else {
      add_to_new_buffer = false;
      if (buffer.mid(i, end.size()) == end) {
        comment = false;
        i += end.size() - 1;
      }
    }
    if (add_to_new_buffer) {
      new_buffer += buffer[i];
    }
    ++i;
  }
  buffer = new_buffer;
}

bool ConfigMap::assign(QString key, QString type, QString value)
{
  if (type == "real") {
    if (typeid(real()) == typeid(float())) {
      type = "float";
    } else {
      type = "double";
    }
  }
  QVariant variant;
  if (type == "float") {
    variant = value.toFloat();
  } else if (type == "double") {
    variant = value.toDouble();
  } else if (type == "int") {
    variant = value.toInt();
  } else if (type == "bool") {
    bool bool_value;
    value = value.toLower();
    if      (value == "true")  bool_value = true;
    else if (value == "on")    bool_value = true;
    else if (value == "yes")   bool_value = true;
    else if (value == "false") bool_value = false;
    else if (value == "off")   bool_value = false;
    else if (value == "no")    bool_value = false;
    else return false;
    variant = bool_value;
  } else if (type == "string") {
    variant = value;
  }
  m_Map[key] = variant;
  return true;
}

void ConfigMap::addFile(QString file_name)
{
  QString buffer;
  QFile file(file_name);
  file.open(QIODevice::ReadOnly);
  buffer = file.readAll();
  processComments(buffer, "/*", "*/");
  processComments(buffer, "//", "\n");
  QStringList items = buffer.split(";");
  foreach (QString item, items) {
    if (item.split("=").size() == 2) {
      QStringList left_side = (item.split("=")[0].trimmed()).split(QRegExp("\\s+"), QString::SkipEmptyParts);
      QString value = item.split("=")[1].trimmed();
      if (left_side.size() != 2) {
        QString msg = "Syntax error : \"" + item + "\"";
        ERROR(qPrintable(msg));
      }
      QString type = left_side[0].trimmed();
      QString key  = left_side[1].trimmed();
      if (!assign(key, type, value)) {
        QString msg = "Syntax error : \"" + item + "\"";
        ERROR(qPrintable(msg));
      }
    }
  }
}

void ConfigMap::addDirectory(QString path)
{
  path = QDir::currentPath() + "/" + path;
  QDir dir(path);
  if (!dir.exists()) {
    QString msg = "The directory \"" + path + "\" does not exists!";
    ERROR(qPrintable(msg));
  }
  QStringList filters;
  filters << "*.dnc";
  filters << "*.Dnc";
  filters << "*.dNc";
  filters << "*.dnC";
  filters << "*.DNc";
  filters << "*.DnC";
  filters << "*.dNC";
  filters << "*.DNC";
  QStringList config_files = dir.entryList(filters, QDir::Files);
  foreach (QString file_name, config_files) {
    addFile(dir.absolutePath() + "/" + file_name);
  }
}

bool ConfigMap::exists(QString key)
{
  if (m_Map.find(key) != m_Map.end()) {
    return true;
  }
  return false;
}
