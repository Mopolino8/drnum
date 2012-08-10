// !!
// This is a part of MOUSE, a library for PDE's on unstructured grids
// Copyright (C) 1999 Oliver Gloth <oliver@vug.uni-duisburg.de>
// Institut fuer Verbrennung und Gasdynamik
// Universitaet Duisburg, Germany
// Institute for Combustion and Gas Dynamics
// University of Duisburg, Germany
// Tue Oct 26 1999
//
// please see http://www.vug.uni-duisburg.de/MOUSE for more information
// please send any questions or suggestions to mouse@www.vug.uni-duisburg.de
//  
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software 
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// !!
#ifndef MError_HH
#define MError_HH
 
#include <namespace_mouse.hh>
BEGIN_MOUSE
 
 
 
END_MOUSE
 
#include <cstddef>
#include <iostream>
 
class MError;
 
#include <mouse_types.hh>
#include <string>
 
BEGIN_MOUSE
 
 
class MError
{
	string name;
 
   string message;
 
 protected:
   virtual void SetName();
 
public:
  MError(string a_message);
 
  virtual ~MError();
 
  string Name();
 
  string Message() { return message; };
};
 
 

END_MOUSE
 
#endif

