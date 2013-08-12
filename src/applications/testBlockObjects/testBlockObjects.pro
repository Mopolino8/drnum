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
# + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
# +                                                                      +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#-------------------------------------------------
#
# Project created by QtCreator 2013-07-27T05:11:42
#
#-------------------------------------------------
TEMPLATE = app
CONFIG += console

drnum_app.path  = ../../../bin
drnum_app.files = testBlockObjects
INSTALLS += drnum_app

include (../drnum_app.pri)

SOURCES      = main.cpp
HEADERS      = main.h \
    geoblockobjecttest001.h \
    geoblockobjecttest002.h \
    geoblockobjecttest003.h \
    geoblockobjecttest004.h \
    geoblockobjecttest005.h \
    geoblockobjecttest006.h

#SOURCES      = main.cpp main.cu
#SOURCES     -= main.cu
#HEADERS      = main.h
#CUDA_SOURCES = main.cu



#QT       += core

#QT       -= gui

#TARGET = testBlockObjects
#CONFIG   += console
#CONFIG   -= app_bundle

#TEMPLATE = app


#SOURCES += main.cpp
