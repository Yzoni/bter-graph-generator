project(bter_seq)

set(HEADER_FILES
        BterSetup.h)

set(SOURCE_FILES
        BterSetup.cpp)

add_library(bter_seq STATIC ${SOURCE_FILES} ${HEADER_FILES})