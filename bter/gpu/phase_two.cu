#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <stdio.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
    if (code != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

__global__ void setup_random_kernel(curandState *state, unsigned long seed) {
    int idx = blockIdx.x * threadIdx.x * blockDim.x;
    curand_init(seed, idx, 0, &state[idx]);
}

__global__ void phase_two_fill(double *phase_two_shift_fill, double *phase_two_sz_fill, curandState *state) {
    int idx = blockIdx.x * threadIdx.x * blockDim.x;
    phase_two_fill[idx] = phase_two_shift_fill[idx] + floor(curand_uniform(&state[idx] * phase_two_sz_fill[idx]);
} 

__global__ void phase_two_bulk(double *phase_two_shift_bulk, double *phase_two_sz_bulk, double *phase_two_bulk, curandState *state) {
    int idx = blockIdx.x * threadIdx.x * blockDim.x;
    phase_two_bulk[idx] = phase_two_shift_bulk[idx] + floor(curand_uniform(&state[idx]) * phase_two_sz_bulk[idx]); 
}

__global__ void phase_two_d(double *phase_two_fill, double *phase_two_bulk) {
    int idx = blockIdx.x * threadIdx.x * blockDim.x;
    if (curand_uniform(&state[idx]) < phase_two_fill[idx]) {
        phase_two[idx] = (int) __double2int_rn(phase_two_fill[idx]) ;
    } else {
        phase_two[idx] = (int) __double2int_rn(phase_two_bulk[idx]); 
    }
}

void cuda_wrapper_phase_two(int phase_two_i, int phase_two_j,
                            int length) {
            
        curandState *devStates;
        cudaMalloc(&devStates, length * sizeof(curandState);        

        int *cuda_i, *cuda_j;

        double *cuda_shift_fill, *cuda_sz_fill; 
        double *cuda_shift_bulk, *cuda_sz_bulk;

        size_t size_output = length * sizeof(int);
        size_t size_input = length * sizeof(double);

        gpuErrchk(cudaMalloc((void **) &cuda_i, size_output));
        gpuErrchk(cudaMalloc((void **) &cuda_j, size_output));

        gpuErrchk(cudaMalloc((void **) &cuda_shift_fill, size_input));
        gpuErrchk(cudaMalloc((void **) &cuda_sz_fill, size_input));

        gpuErrchk(cudaMalloc((void **) &cuda_shift_bulk, size_input));
        gpuErrchk(cudaMalloc((void **) &cuda_sz_bulk, size_input));

        phase_two_fill(phase_two_shift_fill, phase_two_sz_fill, devStates);
        phase_two_bulk(phase_two_shift_bulk, phase_two_sz_bulk, phase_two_bulk, devStates)
       
}
