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
#ifndef SHAREDMEMORY_H
#define SHAREDMEMORY_H

class SharedMemory;

#include <typeinfo>

#include "ipcobject.h"
#include "mutex.h"

class SharedMemory : public IPCObject
{

private:

  int     m_Size;
  char   *m_Buffer;
  size_t  m_MaxNumArrays;
  size_t  m_MaxNameLength;
  size_t  m_ArrayDescrLength;
  size_t  m_Offset;

  template <class T1, class T2> bool TypeMatch() { T1 t1; T2 t2; return typeid(t1) == typeid(t2); }

public:

  enum DataType { Unknown, Real, Integer, Character };

protected:

  int    indexOfName(int i)          { return m_Offset + i * m_ArrayDescrLength; }
  int    indexOfArrayStart(int i)    { return m_Offset + i * m_ArrayDescrLength + m_MaxNameLength; }
  int    indexOfArrayLength(int i)   { return m_Offset + i * m_ArrayDescrLength + m_MaxNameLength + sizeof(int); }
  int    indexOfArrayDataType(int i) { return m_Offset + i * m_ArrayDescrLength + m_MaxNameLength + sizeof(int) + sizeof(DataType); }
  int    arrayStart(int i)           { return *((int*) (&m_Buffer[indexOfArrayStart(i)])); }
  size_t dataIndex();

  template <class T> void get_pointer(int i, T *&t) { t = (T*)(&m_Buffer[i]); }

public:

  SharedMemory(int id_num, int size = 0, bool is_owner = false);
  virtual ~SharedMemory();

  virtual void close();

  int  size()    { return m_Size; }
  void reset();

  int         arrayIndex(std::string name);
  int         arrayLength(int i)          { return *((int*) (&m_Buffer[indexOfArrayLength(i)])); }
  DataType    arrayDataType(int i)        { return *((DataType*) (&m_Buffer[indexOfArrayDataType(i)])); }
  std::string arrayName(int i);
  size_t      numArrays();

  template <class T> void writeArray(std::string name, int length, T *array);
  template <class T> void readArray(std::string name, T *&array);
  template <class T> void writeValue(std::string name, T *value) { writeArray(name, 1, value); }
  template <class T> void readValue(std::string name, T &value);

  std::string readString(std::string name);

};

template <class T>
void SharedMemory::writeArray(std::string name, int length, T *array)
{
  if (name.size() >= m_MaxNameLength) {
    error("SharedMemory::writeArray: 'array name too long'");
  }

  // check if array exists
  int i_array = arrayIndex(name);

  if (i_array >= 0) {

    // compare the data type
    DataType data_type = Unknown;
    if (TypeMatch<T,double>()) data_type = Real;
    if (TypeMatch<T,int>())    data_type = Integer;
    if (TypeMatch<T,char>())   data_type = Character;
    if (data_type != arrayDataType(i_array)) {
      error("SharedMemory::writeArray: 'field \"" + name + "\" exists but has a different type'");
    }
    if (length != arrayLength(i_array)) {
      error("SharedMemory::writeArray: 'field \"" + name + "\" exists but has a different length'");
    }

    // get last data index (start of array)
    T *shm_array = 0;
    get_pointer(arrayStart(i_array), shm_array);

    // write the array
    for (int i = 0; i < length; ++i) {
      shm_array[i] = array[i];
    }

  } else {

    // check if maximal number of arrays is exceeded
    if (numArrays() >= m_MaxNumArrays) {
      error("SharedMemory::writeArray: 'array limit exceeded'");
    }

    // write the name
    for (size_t i = 0; i < m_MaxNameLength; ++i) {
      char c = char(0);
      if (i < name.size()) {
        c = name[i];
      }
      m_Buffer[indexOfName(numArrays()) + i] = c;
    }

    // write the length
    int *L = 0;
    get_pointer(indexOfArrayLength(numArrays()), L);
    *L = length;

    // write the data type
    DataType data_type = Unknown;
    if (TypeMatch<T,double>()) data_type = Real;
    if (TypeMatch<T,int>())    data_type = Integer;
    if (TypeMatch<T,char>())   data_type = Character;
    DataType *DT;
    get_pointer(indexOfArrayDataType(numArrays()), DT);
    *DT = data_type;

    // get and write the data start
    int *DS = 0;
    get_pointer(indexOfArrayStart(numArrays()), DS);
    *DS = dataIndex();

    // get last data index (start of array)
    int *DI = 0;
    get_pointer(4, DI);
    T *shm_array = 0;
    get_pointer(*DI, shm_array);

    // update the last data index
    *DI += length*sizeof(T);
    if (*DI >= m_Size) {
      error("SharedMemory::writeArray: 'buffer size exceeded'");
    }

    // write the array
    for (int i = 0; i < length; ++i) {
      shm_array[i] = array[i];
    }

    // increment the last array index
    int *array_index = 0;
    get_pointer(0, array_index);
    ++(*array_index);
  }
}

template <class T>
void SharedMemory::readArray(std::string name, T *&array)
{
  int i = arrayIndex(name);
  if (i < 0) {
    error("SharedMemory::readArray: 'field \"" + name + "\" not found'");
  }

  // compare the data type
  DataType data_type = Unknown;
  if (TypeMatch<T,double>()) data_type = Real;
  if (TypeMatch<T,int>())    data_type = Integer;
  if (TypeMatch<T,char>())   data_type = Character;
  if (data_type != arrayDataType(i)) {
    error("SharedMemory::readArray: type mismatch for 'field \"" + name + "\"'");
  }


  int L = arrayLength(i);
  array = new T[L];
  T *shm_array = 0;
  get_pointer(arrayStart(i), shm_array);
  for (int j = 0; j < L; ++j) {
    array[j] = shm_array[j];
  }

}

template <class T>
void SharedMemory::readValue(std::string name, T &value)
{
  T *pointer;
  readArray(name, pointer);
  value = *pointer;
  delete [] pointer;
}

#endif // SHAREDMEMORY_H
