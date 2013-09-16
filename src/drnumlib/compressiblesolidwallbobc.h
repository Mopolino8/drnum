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
#ifndef COMPRESSIBLESOLIDWALLBOBC_H
#define COMPRESSIBLESOLIDWALLBOBC_H

class CompressibleSolidWallBOBC;

#include "blockobjectbc.h"

class CompressibleSolidWallBOBC : public BlockObjectBC
{
public:
    CompressibleSolidWallBOBC(size_t field);

    CompressibleSolidWallBOBC(size_t field, BlockObject* block_object);

    virtual void operator()();

};

#endif // COMPRESSIBLESOLIDWALLBOBC_H
