cmake_minimum_required(VERSION 2.8.11)

set(CMAKE_CXX_COMPILER "g++")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11")

project(BTER_sequential)

find_package(CUDA)

set(SOURCE_FILES main.cpp)
add_executable(bter_seq_run ${SOURCE_FILES})

include_directories(bter/gpu bter)

add_subdirectory(bter/gpu)
add_subdirectory(bter)
add_subdirectory(test)
target_link_libraries(bter_seq_run bter)