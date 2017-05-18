cmake_minimum_required(VERSION 3.7)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

project(BTER_sequential)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp)
add_executable(bter_seq_run ${SOURCE_FILES})

include_directories(bter_seq)

add_subdirectory(bter_seq)
add_subdirectory(test)
target_link_libraries(bter_seq_run bter_seq)