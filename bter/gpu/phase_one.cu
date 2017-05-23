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

__global__ void phase_one_shift(double *block_b, double *block_i, double *block_n, int *shift, curandState *state) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    // Compute block and compute its offset
    shift[idx] = (int) __double2int_rn(
            block_i[idx] + __double2int_rd(curand_uniform(&state[idx]) * block_b[idx]) * block_n[idx]);

}

__global__ void phase_one_i(int *i, double *block_b, double *block_i, double *block_n, int *shift, curandState *state) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    // Choose first node
    i[idx] = (int) __double2int_rn(__double2int_rd(curand_uniform(&state[idx]) * block_n[idx]) + shift[idx]);

}

__global__ void phase_one_j(int *i, int *j, double *block_b, double *block_i, double *block_n, int *shift,
                            curandState *state) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Choose second node
    // "Without replacement"
    j[idx] = (int) __double2int_rn(__double2int_rd(curand_uniform(&state[idx]) * (block_n[idx] - 1)) + shift[idx]);

    // Remove loops
    if (j[idx] >= i[idx]) {
        ++j[idx];
    }
}

void cuda_wrapper_phase_one(int *i, int *j,
                            double *block_b, double *block_i, double *block_n,
                            int length) {
    curandState *devStates;
    cudaMalloc(&devStates, length * sizeof(curandState));

    int *cuda_i, *cuda_j;

    double *cuda_block_b, *cuda_block_i, *cuda_block_n;

    size_t size_output = length * sizeof(int);
    size_t size_input = length * sizeof(double);

    gpuErrchk(cudaMalloc((void **) &cuda_i, size_output));
    gpuErrchk(cudaMalloc((void **) &cuda_j, size_output));

    gpuErrchk(cudaMalloc((void **) &cuda_block_b, size_input));
    gpuErrchk(cudaMalloc((void **) &cuda_block_i, size_input));
    gpuErrchk(cudaMalloc((void **) &cuda_j, size_input));

    gpuErrchk(cudaMemcpy(cuda_block_b, block_b, size_input, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(cuda_block_i, block_i, size_input, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(cuda_block_n, block_n, size_input, cudaMemcpyHostToDevice));

    // Temporary shift array only on device
    int *cuda_shift;
    gpuErrchk(cudaMalloc((void **) &cuda_shift, size_input));

    int blocksize = 256;
    int nblock = length / blocksize + (length % blocksize == 0 ? 0 : 1);
    setup_random_kernel <<< nblock, blocksize >>> (devStates, time(NULL));
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    phase_one_shift << < nblock, blocksize >> > (cuda_block_b, cuda_block_i, cuda_block_n, cuda_shift, devStates);

    phase_one_i << < nblock, blocksize >> > (cuda_i, cuda_block_b, cuda_block_i, cuda_block_n, cuda_shift, devStates);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    phase_one_j << < nblock, blocksize >> >
                             (cuda_i, cuda_j, cuda_block_b, cuda_block_i, cuda_block_n, cuda_shift, devStates);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    cudaMemcpy(i, cuda_i, size_output, cudaMemcpyDeviceToHost);
    cudaMemcpy(j, cuda_j, size_output, cudaMemcpyDeviceToHost);

    // FREE
    cudaFree(cuda_shift);

    cudaFree(cuda_i);
    cudaFree(cuda_j);

    cudaFree(cuda_block_b);
    cudaFree(cuda_block_i);
    cudaFree(cuda_block_n);
}

void cuda_wrapper_phase_two() {

}
