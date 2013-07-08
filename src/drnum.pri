CONFIG += cuda

QMAKE_CXXFLAGS += -Wno-deprecated
QMAKE_CXXFLAGS += -fopenmp
QMAKE_CXXFLAGS_RELEASE += -O3
QMAKE_CXXFLAGS_RELEASE += -finline-limit=100000
QMAKE_CXXFLAGS_RELEASE += --param large-function-growth=100000
QMAKE_CXXFLAGS_RELEASE += --param inline-unit-growth=100000

INCLUDEPATH += $(VTKINCDIR)



