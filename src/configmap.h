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
