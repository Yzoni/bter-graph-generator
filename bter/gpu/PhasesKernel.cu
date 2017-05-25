#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <stdio.h>

__global__ void
setup_random_kernel(curandState *state, unsigned long seed, int length) {
    int idx = blockIdx.x * threadIdx.x * blockDim.x;
    if (idx < length) {
        curand_init(seed, idx, 0, &state[idx]);
    }
}

/*
 *  PHASE ONE
 */

__global__ void
phase_one_shift(double *block_b, double *block_i, double *block_n, int *shift, curandState *state, int length) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    // Compute block and compute its offset
    if (idx < length) {
        shift[idx] = (int) __double2int_rn(
                block_i[idx] + __double2int_rd(curand_uniform(&state[idx]) * block_b[idx]) * block_n[idx]);
    }
}

__global__ void
phase_one_i(int *i, double *block_b, double *block_i, double *block_n, int *shift, curandState *state, int length) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    // Choose first node
    if (idx < length) {
        i[idx] = (int) __double2int_rn(__double2int_rd(curand_uniform(&state[idx]) * block_n[idx]) + shift[idx]);
    }
}

__global__ void
phase_one_j(int *i, int *j, double *block_b, double *block_i, double *block_n, int *shift,
            curandState *state, int length) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < length) {
        // Choose second node
        // "Without replacement"
        j[idx] = (int) __double2int_rn(__double2int_rd(curand_uniform(&state[idx]) * (block_n[idx] - 1)) + shift[idx]);

        // Remove loops
        if (j[idx] >= i[idx]) {
            ++j[idx];
        }
    }
}

/*
 *  PHASE TWO
 */

__global__ void
phase_two_fill(double *phase_two_shift_fill, double *phase_two_sz_fill, double *phase_two_fill,
               curandState *state, int length) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < length) {
        phase_two_fill[idx] =
                phase_two_shift_fill[idx] + __double2int_rd(curand_uniform(&state[idx]) * phase_two_sz_fill[idx]);
    }
}

__global__ void
phase_two_bulk(double *phase_two_shift_bulk, double *phase_two_sz_bulk, double *phase_two_bulk,
               curandState *state, int length) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < length) {
        phase_two_bulk[idx] =
                phase_two_shift_bulk[idx] + __double2int_rd(curand_uniform(&state[idx]) * phase_two_sz_bulk[idx]);
    }
}

__global__ void
phase_two_d(double *phase_two_fill, double *phase_two_bulk, int *phase_two,
            curandState *state, int length) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < length) {
        if (curand_uniform(&state[idx]) < phase_two_fill[idx]) {
            phase_two[idx] = (int) __double2int_rn(phase_two_fill[idx]);
        } else {
            phase_two[idx] = (int) __double2int_rn(phase_two_bulk[idx]);
        }
    }
}