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
#ifndef CODESTRING_H
#define CODESTRING_H

//#include <cstddef>
//#include <string.h>
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <vector>

#include "drnum.h"

#include <iostream>
#include <sstream>
#include <string>
using namespace std;

//class CodeString;

class CodeString : public string
{

  // string m_mystring;    ///< String containing code words. Code words: any ascii string

protected: // methods

public: // methods

  /**
    * Empty constructor.
    */
  CodeString();

  /**
   * Constructor upon a given string
   * @param inputstring the string to use
   */
  CodeString(const string& inputstring);

  /**
   * Constructor upon an input stream
   * @param inputstring the string to use
   */
  CodeString(istringstream& input_iss);

  /**
    * Build up.
    * inputiss input stream
    */
  void buildFrom(istringstream& inputiss);

  /**
    * Comparison method.
    * Note: concatcodestring is invalid, if codestring do not match.
    * @param otherstring string or codestring to compare to
    * @param allsame true, if both match identically (return reference)
    * @param overruled true, if one or more code words of "this" are "0",
    *                  while corresponding in otherstring arent (ret ref).
    * @param underruled true, if one or more code words of otherstring
    *                   are "0", while corresponding in "this" arent (ret ref).
    * @return concatcodestring concatenated codestring replacing all underruled and overruled
    *                          "0"es by corresponding non-zero code words of this or otherstring
    */
  CodeString compareCodeString(const CodeString& otherstring,
                               bool& allsame,
                               bool& overruled,
                               bool& underruled);


  /** Count
    * @return number of code words in string.
    */
  size_t countNumCodeWords();


  /** Access a code word
    * @param i_word the i_th code word to extract.
    * @param success true, if i_th word exists (return reference)
    * @return the code word as a string.
    */
  string accessCodeWord(const size_t& i_word,
                        bool& success);


  /**
    * Destructor. No data beyond string though.
    */
  virtual ~CodeString();

};

#endif // CODESTRING_H
