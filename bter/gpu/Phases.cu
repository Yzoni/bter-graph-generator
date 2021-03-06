#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <stdio.h>

#include "PhasesKernel.cu"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
    if (code != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

void cuda_setup_random_kernel(int length, curandState *devStates) {
    int blocksize = 512;
    int nblock = length / blocksize + (length % blocksize == 0 ? 0 : 1);

    printf("CuRand init: %d nblock: %d and blocksize: %d\n", length, nblock, blocksize);
    setup_random_kernel << < nblock, blocksize >> > (devStates, length, 0);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());
}

void cuda_wrapper_rand_array(int length, double *out_array) {
    curandState *devStates;
    gpuErrchk(cudaMalloc(&devStates, length * sizeof(curandState)));
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    double *cuda_rand_array;
    gpuErrchk(cudaMalloc((void **) &cuda_rand_array, length * sizeof(double)));
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    int blocksize = 512;
    int nblock = length / blocksize + (length % blocksize == 0 ? 0 : 1);

    cuda_setup_random_kernel(length, devStates);
    printf("CuRand array: length: %d nblock: %d and blocksize: %d\n", length, nblock, blocksize);

    get_random_array << < nblock, blocksize >> > (devStates, length, cuda_rand_array);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    cudaMemcpy(out_array, cuda_rand_array, length * sizeof(double), cudaMemcpyDeviceToHost);
    printf("CuRand produced (f10): ");
    for (int i = 0; i < 10; ++i) printf("%f ", out_array[i]);
    printf("\n");

    cudaFree(devStates);
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
    gpuErrchk(cudaMalloc((void **) &cuda_block_n, size_input));

    gpuErrchk(cudaMemcpy(cuda_block_b, block_b, size_input, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(cuda_block_i, block_i, size_input, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(cuda_block_n, block_n, size_input, cudaMemcpyHostToDevice));

    // Temporary shift array only on device
    int *cuda_shift;
    gpuErrchk(cudaMalloc((void **) &cuda_shift, size_input));

    int blocksize = 512;
    int nblock = length / blocksize + (length % blocksize == 0 ? 0 : 1);
    cuda_setup_random_kernel(length, devStates);

    // Shift
    phase_one_shift << < nblock, blocksize >> >
                                 (cuda_block_b, cuda_block_i, cuda_block_n, cuda_shift, devStates, length);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    // i
    phase_one_i << < nblock, blocksize >> >
                             (cuda_i, cuda_block_b, cuda_block_i, cuda_block_n, cuda_shift, devStates, length);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    // j
    phase_one_j << < nblock, blocksize >> >
                             (cuda_i, cuda_j, cuda_block_b, cuda_block_i, cuda_block_n, cuda_shift, devStates, length);
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

void cuda_wrapper_phase_two(int *phase_two,
                            double *phase_two_shift_fill, double *phase_two_sz_fill,
                            double *phase_two_shift_bulk, double *phase_two_sz_bulk,
                            double *phase_two_rd_fill,
                            int length) {

    curandState *devStates;
    cudaMalloc(&devStates, length * sizeof(curandState));

    int *cuda_phase_two;

    double *cuda_shift_fill, *cuda_sz_fill;
    double *cuda_shift_bulk, *cuda_sz_bulk;
    double *cuda_rd_fill;

    size_t size_output = length * sizeof(int);
    size_t size_input = length * sizeof(double);

    gpuErrchk(cudaMalloc((void **) &cuda_phase_two, size_output));

    gpuErrchk(cudaMalloc((void **) &cuda_shift_fill, size_input));
    gpuErrchk(cudaMalloc((void **) &cuda_sz_fill, size_input));

    gpuErrchk(cudaMalloc((void **) &cuda_shift_bulk, size_input));
    gpuErrchk(cudaMalloc((void **) &cuda_sz_bulk, size_input));

    gpuErrchk(cudaMalloc((void **) &cuda_rd_fill, size_input));

    double *cuda_fill, *cuda_bulk;
    gpuErrchk(cudaMalloc((void **) &cuda_fill, size_input));
    gpuErrchk(cudaMalloc((void **) &cuda_bulk, size_input));

    gpuErrchk(cudaMemcpy(cuda_shift_fill, phase_two_shift_fill, size_input, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(cuda_sz_fill, phase_two_sz_fill, size_input, cudaMemcpyHostToDevice));

    gpuErrchk(cudaMemcpy(cuda_shift_bulk, phase_two_shift_bulk, size_input, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(cuda_sz_bulk, phase_two_sz_bulk, size_input, cudaMemcpyHostToDevice));

    gpuErrchk(cudaMemcpy(cuda_rd_fill, phase_two_rd_fill, size_input, cudaMemcpyHostToDevice));

    int blocksize = 512;
    int nblock = length / blocksize + (length % blocksize == 0 ? 0 : 1);
    cuda_setup_random_kernel(length, devStates);

    phase_two_fill << < nblock, blocksize >> > (cuda_shift_fill, cuda_sz_fill, cuda_fill, devStates, length);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    phase_two_bulk << < nblock, blocksize >> > (cuda_shift_bulk, cuda_sz_bulk, cuda_bulk, devStates, length);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    phase_two_d << < nblock, blocksize >> > (cuda_fill, cuda_bulk, cuda_phase_two, cuda_rd_fill, devStates, length);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    // COPY BACK
    gpuErrchk(cudaMemcpy(phase_two, cuda_phase_two, size_output, cudaMemcpyDeviceToHost));

    // FREE
    cudaFree(cuda_shift_fill);
    cudaFree(cuda_sz_fill);
    cudaFree(cuda_shift_bulk);
    cudaFree(cuda_sz_bulk);

    cudaFree(cuda_fill);
    cudaFree(cuda_bulk);

    cudaFree(cuda_phase_two);
}
