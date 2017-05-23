project(bter_lib)

set(HEADER_FILES
        BterSetup.h
        BterPhases.h
        )

set(SOURCE_FILES
        BterSetup.cpp
        BterPhases.cpp
        )

add_library(bter STATIC ${SOURCE_FILES} ${HEADER_FILES})

target_link_libraries(bter bter_gpu)
