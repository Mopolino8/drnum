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
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef STRINGTOOLS_HH
#define STRINGTOOLS_HH

#include <string>
#include <vector>
#include <list>
#include <sstream>
#include <iostream>
#include <algorithm>

using namespace std;

namespace StringTools
{


/**
 * append a character to the left side of a string until a given length has been reached.
 * @param s the source string
 * @param c the character to append
 * @param l the desired length
 * @return the filled string
 */
inline string leftFill(string s, char c, size_t l)
{
  while (s.size() < l) {
    s = c + s;
  }
  return s;
}


/**
 * append a character to the right side of a string until a given length has been reached.
 * @param s the source string
 * @param c the character to append
 * @param l the desired length
 * @return the filled string
 */
inline string rightFill(string s, char c, size_t l)
{
  while (s.size() < l) {
    s = s + c;
  }
  return s;
}


/**
 * convert a string into another type.
 * <b>Attention</b> this method only works for a type <i>T</i> for
 * which an input operator <i>operator>></i>exists.
 * @param s the string to convert
 * @param t a variable which will contain the converted value
 * @return true if the operation succeeded
 */
template <typename T>
bool stringTo(string s, T &t)
{
  s += " -";
  istringstream stream(s);
  stream >> t;
  return stream.good();
}


/**
 * convert another type into a string.
 * <b>Attention</b> this method only works for a type <i>T</i> for
 * which an output operator <i>operator<<</i>exists.
 * @param t the variable to convert
 * @param s the resulting string
 */
template <typename T>
void toString(T t, string &s)
{
  ostringstream stream;
  stream << t;
  s = stream.str();
}


/**
 * convert another type into a string.
 * <b>Attention</b> this method only works for a type <i>T</i> for
 * which an output operator <i>operator<<</i>exists.
 * @param t the variable to convert
 * @return the resulting string
 */
template <typename T>
string toString(T t, int fill = 0)
{
  ostringstream stream;
  stream << t;
  string s = stream.str();
  if      (fill < 0) s = leftFill(s,' ',-fill);
  else if (fill > 0) s = rightFill(s,' ',fill);
  return s;
}


/**
 * replace a character with a different one.
 * @param s the source string
 * @param c_orig the character to replace
 * @param c_new the new character
 * return a string with the replaced characters
 */
inline string replace(string s, char c_orig, char c_new)
{
  string s_new = "";
  string::iterator i = s.begin();
  while (i != s.end()) {
    if (*i == c_orig) s_new += c_new;
    else s_new += *i;
    i++;
  }
  return s_new;
}


/**
 * read an entire line from an input stream
 * @param s the stream to read from
 * @return the line read
 */
inline string readLine(istream &s)
{
  string line = "";
  bool done = false;
  do {
    char c;
    s.get(c);
    if (c == '\n') {
      done = true;
    } else {
      line += c;
      done = (s.eof() || s.fail());
    }
  } while (!done);
  return line;
}


/**
 * extract a part of a string.
 * @param s the original string
 * @param i1 index of the first character of the intended sub-string
 * @param i2 index of the last character of the intended sub-string
 * @return the sub-string
 */
inline string subString(string s, size_t i1, size_t i2)
{
  string sub = "";
  size_t i = i1;
  while ((i <= i2) && (i < s.size())) {
    sub += s[i];
    ++i;
  };
  return sub;
}


/**
 * extract a part at the end of a string.
 * @param s the original string
 * @param n number of characters of the intended sub-string
 * @return the sub-string
 */
inline string right(string s, size_t n)
{
  return subString(s,s.size()-n-1,s.size()-1);
}


/**
 * extract a part at the beginning of a string.
 * @param s the original string
 * @param n number of characters of the intended sub-string
 * @return the sub-string
 */
inline string left(string s, size_t n)
{
  return subString(s,0,n-1);
}


/**
 * split a string into several strings.
 * @param s the string to be split
 * @param sub the list which the sub strings will be appeneded to
 * @param c the dividing character
 * @return the number of sub strings
 */
inline int split(string s, list<string> &sub, char c = ' ')
{
  string word = "";
  bool first = true;
  int N = 0;
  for (size_t i = 0; i < s.size(); ++i) {
    if (s[i] != c) {
      first = false;
      word += s[i];
    } else {
      if (!first) {
        sub.push_back(word);
        ++N;
        word = "";
        first = true;
      }
    }
  }
  if (word.size() > 0) {
    sub.push_back(word);
    ++N;
  }
  return N;
}


/**
 * split a string into several strings.
 * @param s the string to be split
 * @param sub the vector which will contain the sub strings
 * @param c the dividing character
 * @return the number of sub strings
 */
inline int split(string s, vector<string> &sub, char c = ' ')
{
  list<string> l;
  int N = split(s,l,c);
  sub.resize(N);
  copy(l.begin(),l.end(),sub.begin());
  return N;
}


}

#endif


