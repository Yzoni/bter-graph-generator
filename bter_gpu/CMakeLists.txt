project(bter_gpu)

cuda_add_library(bter_gpu test.cu)

target_link_libraries(bter_gpu ${CUDA_LIBRARIES})