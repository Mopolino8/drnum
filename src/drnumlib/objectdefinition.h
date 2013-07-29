#ifndef OBJECTDEFINITION_H
#define OBJECTDEFINITION_H

class ObjectDefinition;

#include "blockcfd.h"

class ObjectDefinition
{
  /** @todo There seams to be a bug in the C++ compiler. It will not
   *  accept other instances of derived classes of ObjectDefinition
   *  access the below functions, if in protected range (???) */

  //C++ bug? friend class ObjectDefinition;

protected: // attributes
  bool m_KnownInside;

protected: // methods
  //C++ bug? setKnownInside(const bool& inside) {m_KnownInside = inside;}

  /** @todo define wether to use o-coords, or generalized parent coords xyz
    *       in all derived class. */

public:

  void setKnownInside(const bool& inside) {m_KnownInside = inside;}
  bool getKnownInside() {return m_KnownInside;}

  ObjectDefinition();

  virtual bool isInside(const real& xo, const real& yo, const real& zo) = 0;

  bool isInside(vec3_t xyzo) { return isInside (xyzo[0], xyzo[1],xyzo[2]); }

  virtual void getLowestObjects(vector<ObjectDefinition*>& my_lowest_objects);

  virtual bool evalBool() {return getKnownInside();}

};

#endif // OBJECTDEFINITION_H
