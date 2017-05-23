project(bter_gpu)

cuda_add_library(bter_gpu phase_one.cu phase_two.cu)

target_link_libraries(bter_gpu ${CUDA_LIBRARIES})
