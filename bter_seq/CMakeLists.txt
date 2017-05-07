project(bter_seq)

set(HEADER_FILES
        BterSetup.h
        BterPhases.h)

set(SOURCE_FILES
        BterSetup.cpp
        BterPhases.cpp)

add_library(bter_seq STATIC ${SOURCE_FILES} ${HEADER_FILES})