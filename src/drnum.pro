TEMPLATE  = subdirs
LANGUAGE  = C++
CONFIG   += ordered recursive

include (drnum.pri)

SUBDIRS  += drnumlib
SUBDIRS  += applications

drnumlib.file    = drnumlib/drnumlib.pro

applications.file    = applications/applications.pro
applications.depends = drnumlib
