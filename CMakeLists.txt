cmake_minimum_required(VERSION 2.8.11)

set(CMAKE_CXX_COMPILER "g++")
add_definitions(-std=c++0x)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++0x")

# 5.2 compatibility for Titan X
# 6.1 compatibility for GTX 1070
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11
--generate-code arch=compute_52,code=sm_52
--generate-code arch=compute_61,code=sm_61")


project(bter_project)

find_package(CUDA)

set(PYTHON_VERSION 3.4)
set(PYTHON_INCLUDE_DIR "/cm/shared/package/python/3.5.2/include/python3.5m")
set(PYTHON_LIBRARY "/cm/shared/package/python/3.5.2/lib/libpython3.5m.so.1.0")

set(Python_ADDITIONAL_VERSIONS 3.4)
find_package(PythonLibs 3.4 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

set(SOURCE_FILES main.cpp)
add_executable(bter_run ${SOURCE_FILES})

include_directories(bter bter/gpu spdlog)
link_directories(spdlog)

add_subdirectory(bter/gpu)
add_subdirectory(bter)
add_subdirectory(test)

target_link_libraries(bter_run bter ${PYTHON_LIBRARIES})

file(COPY parameters DESTINATION ${CMAKE_BINARY_DIR})