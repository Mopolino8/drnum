// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
// +                                                                      +
// + enGrid is free software: you can redistribute it and/or modify       +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 
#ifndef UTILITIES_H
#define UTILITIES_H

/** @file utilities.h Contains several utility functions. */

#include <QVector>
#include <QMap>
#include <QObject>
#include <QtDebug>

#include <iostream>
using namespace std;

#include "math/mathvector.h"
#include "math/smallsquarematrix.h"
#include "math/linsolve.h"

using namespace std;

/** Restricts value to the CYCLIC [min,max[ range.
 * Equivalent to min + modulo(value-min, max-min)
 * @param value value to restrict
 * @param min minimum
 * @param max maximum (if value=max, then min will be returned)
 */
inline real restrict_to(real value, real min, real max)
{
  return value - floor((value - min) / (max - min))*(max - min);
}

/// Converts degrees to radians
inline real RadiansFromDegrees(real a_deg)
{
  return(a_deg / 180.0*M_PI);
}

/// Converts radians to degrees
inline real DegreesFromRadians(real a_rad)
{
  return(a_rad / M_PI*180.0);
}

real toreal(QString str);///< Equivalent of the QString::toreal function, but it also replaces "," with "."

/** Converts a real to a string.
 * @param x real value to convert
 * @param separator symbol to use as decimal separator. It is recommended to always use the default: ".".
 */
QString toString(real x, QString separator = QObject::tr("."));

/// a modulo function which also works for negative numbers
inline int modulo(int a, int b)
{
  return((a % b + b) % b);
}

/// a modulo function which also works for negative numbers
inline real modulo(real a, real b)
{
  while (a < 0) a += b;
  while (a >= b) a -= b;
  return a;
}

/** Adds a suffix to a filename if it is missing.
 * @param filename string to which to add the suffix
 * @param suffix suffix to add. suffix should be of the form 'ext' without the dot!!!
 * @param remove_old_suffix If true, any previous suffix of filename will be removed. ex: foo.bar -> foo.png
 * @return the filename with the new suffix
 */
QString addSuffix(QString filename, QString suffix, bool remove_old_suffix);

///////////////////////////////////////////

/** Transposes a vector of vectors (matrix)
 * @param in the matrix to transpose
 * @param nrows the number of rows of the input matrix
 * @param ncolumns the number of columns of the input matrix
 * @return the transposed "matrix"
 */
template <class T>
QVector < QVector <T> > transpose(QVector < QVector <T> > & in, int nrows, int ncolumns)
{
  QVector < QVector <T> > out(ncolumns);
  QVector <T> col(nrows);
  out.fill(col);

  for (int i = 0; i < in.size(); i++) {
    QVector <T> row = in[i];
    for (int j = 0; j < row.size(); j++) {
      out[j][i] = in[i][j];
    }
  }
  return(out);
}

///////////////////////////////////////////

/// ostream operator for QVectors
template <class T>
ostream &operator<<(ostream &out, QVector<T> const & vector)
{
  int N = vector.size();
  out << "[";
  for (int i = 0; i < N; ++i) {
    out << vector.at(i);
    if (i != N - 1) out << ",";
  }
  out << "]";
  return(out);
}

/// ostream operator for QSets
template <class T>
ostream &operator<<(ostream &out, QSet<T> const & set)
{
  out << "[ ";
  foreach(T value, set) out << value << " ";
  out << "]";
  return(out);
}

/// ostream operator for QVectors of QSets
template <class T>
ostream &operator<<(ostream &out, QVector<QSet<T> > & vector)
{
  int N = vector.size();
  out << "[";
  for (int i = 0; i < N; ++i) {
    QSet<T> set = vector.at(i);
    out << set;
    if (i != N - 1) out << ",";
  }
  out << "]";
  return(out);
}

/// ostream operator for QVectors of QVectors
template <class T>
ostream &operator<<(ostream &out, QVector<QVector<T> > & vector)
{
  int N = vector.size();
  out << "[";
  for (int i = 0; i < N; ++i) {
    QVector<T> subvector = vector.at(i);
    out << subvector;
    if (i != N - 1) out << ",";
  }
  out << "]";
  return(out);
}

/// ostream operator for QMaps
template <class T1, class T2>
ostream &operator<<(ostream &out, QMap<T1, T2> & map)
{
  QMapIterator<T1, T2> i(map);
  out << "[";
  while (i.hasNext()) {
    i.next();
    out << " [" << i.key() << ": " << i.value() << "]";
  }
  out << "]";
  return(out);
}

/// ostream operator for a QVector of pairs
template <class T1, class T2>
ostream &operator<<(ostream &out, QVector < pair<T1, T2> > & vector)
{
  int N = vector.size();
  out << "[";
  for (int i = 0; i < N; ++i) {
    out << "<";
    out << vector.at(i).first;
    out << ",";
    out << vector.at(i).second;
    out << ">";
    if (i != N - 1) out << ",";
  }
  out << "]";
  return(out);
}

////////////////////////////////////////////////////

template <class T>
QVector <T> set2Vector(QSet <T> a_set, bool a_sort)
{
  QVector <T> l_vector(a_set.size());
  qCopy(a_set.begin(),a_set.end(),l_vector.begin());
  if(a_sort) qSort(l_vector.begin(),l_vector.end());
  return(l_vector);
}

template <class T>
QSet <T> vector2Set(QVector <T> a_vector, bool a_sort)
{
  QSet <T> l_set;
  foreach(T element, a_vector) l_set.insert(element);
  if(a_sort) qSort(l_set.begin(),l_set.end());
  return(l_set);
}

template <class T>
bool duplicates(QVector <T> a_vector)
{
  QSet <T> l_set;
  foreach(T element, a_vector) l_set.insert(element);
  return l_set.size()!=a_vector.size();
}

/// V(rotated base) = rotationMatrix_X * V(original base)
mat3_t rotationMatrixX(real a_rad);

/// V(rotated base) = rotationMatrix_Y * V(original base)
mat3_t rotationMatrixY(real a_rad);

/// V(rotated base) = rotationMatrix_Z * V(original base)
mat3_t rotationMatrixZ(real a_rad);

/// V(rotated base) = rotateRelativeZXZ * V(original base)
mat3_t rotateRelativeZXZ(real angle_1_rad, real angle_2_rad, real angle_3_rad);

/// V(rotated base) = rotateAbsoluteZXY * V(original base)
mat3_t rotateAbsoluteZXY(real angle_1_rad, real angle_2_rad, real angle_3_rad);

/// returns the polar angle in radians
real getGamma(vec3_t V);

/// returns the  azimuth angle in radians. returns phi from -pi to pi
real getPhi(vec3_t V);

/// A version of QFileDialog::getExistingDirectory which allows preselecting a directory
QString getDirectory(QWidget * parent = 0, const QString & caption = QString(), const QString & selected = QString());
// QString getDirectory( QWidget * parent = 0, const QString & caption = QString(), const QString & dir = QString(), const QString & selected = QString() , Options options = ShowDirsOnly );

bool checkVector(vec3_t V);
bool checkVector(vec2_t V);

// x^3 + a x^2 + b x + c = 0
int polySolveCubic(real a, real b, real c, real * x0, real * x1, real * x2);

// a x^2 + b x + c = 0
int polySolveQuadratic(real a, real b, real c, real * x0, real * x1);

#endif
