#include "configmap.h"

#include <iostream>
#include <fstream>

ConfigMap::ConfigMap()
{
}

void ConfigMap::readFromFile(string file_name)
{
  using namespace std;

  ifstream f(file_name.c_str());
  string buffer;
  string line;
  while (f){
    getline(f, line);
    buffer += line;
  }
  cout << buffer << endl;
}
