TEMPLATE = app
CONFIG += console

drnum_app.path  = ../../../bin
drnum_app.files = drnumBasicAero
INSTALLS += drnum_app

include (../drnum_app.pri)

SOURCES      = main.cpp main.cu
SOURCES     -= main.cu
HEADERS      = main.h
CUDA_SOURCES = main.cu


