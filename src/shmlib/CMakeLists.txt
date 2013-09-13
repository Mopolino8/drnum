SET(shmlib_SOURCES
    barrier.cpp
    ipcexception.cpp
    ipcobject.cpp
    mutex.cpp
    semaphore.cpp
    sharedmemory.cpp
    shmlib.cpp
)

ADD_LIBRARY(shmlib SHARED ${shmlib_SOURCES})

SET_TARGET_PROPERTIES(shmlib PROPERTIES VERSION ${DRNUM_VERSION})

INSTALL(TARGETS shmlib LIBRARY DESTINATION lib)
