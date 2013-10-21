# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of DrNUM.                                          +
# +                                                                      +
# + Copyright 2013 numrax GmbH, enGits GmbH                              +
# +                                                                      +
# + DrNUM is free software: you can redistribute it and/or modify        +
# + it under the terms of the GNU General Public License as published by +
# + the Free Software Foundation, either version 3 of the License, or    +
# + (at your option) any later version.                                  +
# +                                                                      +
# + DrNUM is distributed in the hope that it will be useful,             +
# + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
# + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
# + GNU General Public License for more details.                         +
# +                                                                      +
# + You should have received a copy of the GNU General Public License    +
# + along with DrNUM. If not, see <http://www.gnu.org/licenses/>.        +
# +                                                                      +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TEMPLATE     = lib
shmlib.path  = ../../lib
shmlib.files = libshmlib.*
INSTALLS    += shmlib
QT          -= gui core

include (../drnum.pri)

SOURCES += shmlib.cpp
SOURCES += sharedmemory.cpp
SOURCES += semaphore.cpp
SOURCES += ipcobject.cpp
SOURCES += mutex.cpp
SOURCES += barrier.cpp
SOURCES += ipcexception.cpp

HEADERS += sharedmemory.h
HEADERS += semaphore.h
HEADERS += ipcobject.h
HEADERS += mutex.h
HEADERS += barrier.h
HEADERS += sharedmemory.f.h
HEADERS += mutex.f.h
HEADERS += barrier.f.h
HEADERS += ipcexception.h

OTHER_FILES += test_write.f
OTHER_FILES += shock_tube.f
