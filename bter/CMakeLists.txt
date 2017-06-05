project(bter_lib)

set(HEADER_FILES
        BterSetup.h
        BterPhasesSeq.h
        BterPhasesGpu.h
        BterSamples.h
        )

set(SOURCE_FILES
        BterSetup.cpp
        BterPhasesSeq.cpp
        BterPhasesGpu.cpp
        )

add_library(bter STATIC ${SOURCE_FILES} ${HEADER_FILES})

target_link_libraries(bter bter_gpu)
