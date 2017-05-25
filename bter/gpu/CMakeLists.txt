project(bter_gpu)

cuda_add_library(bter_gpu PhasesKernel.cu Phases.cu)

target_link_libraries(bter_gpu ${CUDA_LIBRARIES})
