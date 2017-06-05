#ifndef BTER_BTERPHASESGPU_H
#define BTER_BTERPHASESGPU_H

#include "BterSetup.h"
#include "BterSamples.h"

void cuda_wrapper_rand_array(int length, double *out_array);

void cuda_wrapper_phase_one(int *phase_one_i, int *phase_one_j,
                            double *block_b, double *block_i, double *block_n,
                            int length);

void cuda_wrapper_phase_two(int *phase_two,
                            double *phase_two_shift_fill, double *phase_two_sz_fill,
                            double *phase_two_shift_bulk, double *phase_two_sz_bulk,
                            int length);

namespace bter {

    class BterPhasesGpu {
    public:
        BterPhasesGpu(BTERSetupResult *bterSetupResult, int dmax, double *nd, double *cd);

        void computeSamples();

        void phaseOnePrepare(int *group_sample, double *block_b, double *block_i, double *block_n);

        void phaseOne(int *phase_one_i, int *phase_one_j);

        void phaseTwoNodePrepare(double *id_bulk, double *nd_bulk,
                                 double *phase_two_shift_fill, double *phase_two_sz_fill,
                                 double *phase_two_shift_bulk, double *phase_two_sz_bulk);

        void phaseTwo(int *phase_two_i, int *phase_one_j);

//        void removeLoopsPhaseTwo(int *i, int *j, int length, std::vector<int, std::allocator<int>> *new_i,
//                                 vector<int, std::allocator<int>> *new_j);

        BterSamples bterSamples;

    private:
        BTERSetupResult *bterSetupResult;

        int dmax;

        double *nd, *cd;

        void randomSample(double *wg, int nsamples, int *binindices);

        void phaseTwoNode(double *id_bulk, double *nd_bulk, int *phase_two);
    };
}

#endif //BTER_BTERPHASESGPU_H
