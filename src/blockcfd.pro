TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    patch.cpp \
    cartesianpatch.cpp \
    compressiblecartesianpatch.cpp

HEADERS += \
    patch.h \
    cartesianpatch.h \
    blockcfd.h \
    compressiblecartesianpatch.h

