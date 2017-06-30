#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <stdio.h>
#include <ctime>

__global__ void
setup_random_kernel(curandState *state, int length, int offset) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < length) {
        curand_init((unsigned long long) clock(), idx, 0, &state[idx]);
    }
}

__global__ void
get_random_array(curandState *state, int length, double *out_array) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < length) {
        curandState localState = state[idx];
        out_array[idx] = curand_uniform_double(&localState);
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
        curandState localState = state[idx];
        shift[idx] = (int) __double2int_rn(
                block_i[idx] + __double2int_rd(curand_uniform(&localState) * block_b[idx]) * block_n[idx]);
        state[idx] = localState;
    }
}

__global__ void
phase_one_i(int *i, double *block_b, double *block_i, double *block_n, int *shift, curandState *state, int length) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    // Choose first node
    if (idx < length) {
        curandState localState = state[idx];
        i[idx] = (int) __double2int_rn(__double2int_rd(curand_uniform(&localState) * block_n[idx]) + shift[idx]);
        state[idx] = localState;
    }
}

__global__ void
phase_one_j(int *i, int *j, double *block_b, double *block_i, double *block_n, int *shift,
            curandState *state, int length) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < length) {
        curandState localState = state[idx];

        // Choose second node
        // "Without replacement"
        j[idx] = (int) __double2int_rn(__double2int_rd(curand_uniform(&localState) * (block_n[idx] - 1)) + shift[idx]);

        // Remove loops
        if (j[idx] >= i[idx]) {
            ++j[idx];
        }
        state[idx] = localState;
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
        curandState localState = state[idx];
        phase_two_fill[idx] =
                phase_two_shift_fill[idx] + __double2int_rd(curand_uniform(&localState) * phase_two_sz_fill[idx]);
        state[idx] = localState;
    }
}

__global__ void
phase_two_bulk(double *phase_two_shift_bulk, double *phase_two_sz_bulk, double *phase_two_bulk,
               curandState *state, int length) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < length) {
        curandState localState = state[idx];
        phase_two_bulk[idx] =
                phase_two_shift_bulk[idx] + __double2int_rd(curand_uniform(&localState) * phase_two_sz_bulk[idx]);
        state[idx] = localState;
    }
}

__global__ void
phase_two_d(double *phase_two_fill, double *phase_two_bulk, int *phase_two, double *phase_two_rd_fill,
            curandState *state, int length) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < length) {
        curandState localState = state[idx];
        if (curand_uniform(&localState) < phase_two_rd_fill[idx]) {
            phase_two[idx] = (int) __double2int_rn(phase_two_fill[idx]);
        } else {
            phase_two[idx] = (int) __double2int_rn(phase_two_bulk[idx]);
        }
        state[idx] = localState;
    }
}