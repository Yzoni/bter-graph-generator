#ifndef BTER_SEQUENTIAL_BTERPHASES_H
#define BTER_SEQUENTIAL_BTERPHASES_H

#include "BterSetup.h"

void cuda_wrapper_rand_array(int length, double *out_array);

void cuda_wrapper_phase_one(int *phase_one_i, int *phase_one_j,
                            double *block_b, double *block_i, double *block_n,
                            int length);

void cuda_wrapper_phase_two(int *phase_two,
                            double *phase_two_shift_fill, double *phase_two_sz_fill,
                            double *phase_two_shift_bulk, double *phase_two_sz_bulk,
                            int length);

namespace bter {

    struct BterSamples {
        int s1, s2;

        BterSamples(int s1, int s2) : s1(s1), s2(s2) {}
    };

    struct BterEdges {
        int *phase_one_i;
        int *phase_one_j;
        int *phase_two_i;
        int *phase_two_j;

        int bter_phase_one_size;
        int bter_phase_two_size;
    };

    class BterPhases {
    public:
        BterPhases(BTERSetupResult *bterSetupResult, int dmax, double *nd, double *cd);

        void computeSamplesSeq();

        void computeSamplesGpu();

        void phaseOnePrepare(int *group_sample, double *block_b, double *block_i, double *block_n);

        void phaseOneSeq(int *phase_one_i, int *phase_one_j);

        void phaseOneGpu(int *phase_one_i, int *phase_one_j);

        void phaseTwoNodePrepare(double *id_bulk, double *nd_bulk,
                                 double *phase_two_shift_fill, double *phase_two_sz_fill,
                                 double *phase_two_shift_bulk, double *phase_two_sz_bulk);

        void phaseTwoSeq(int *phase_two_i, int *phase_two_j);

        void phaseTwoGpu(int *phase_two_i, int *phase_one_j);

        void phaseTwoNodeSeq(double *id_bulk, double *nd_bulk, int *phase_two);

        void phaseTwoNodeGpu(double *id_bulk, double *nd_bulk, int *phase_two);

//        void removeLoopsPhaseTwo(int *i, int *j, int length, std::vector<int, std::allocator<int>> *new_i,
//                                 vector<int, std::allocator<int>> *new_j);

        BterSamples bterSamples;

    private:
        BTERSetupResult *bterSetupResult;

        int dmax;

        double *nd, *cd;

        double randomUnified(int from, int to);

        void randomSample(double *wg, int nsamples, int *binindices);

        void phaseTwoNode(double *id_bulk, double *nd_bulk, int *phase_two);
    };
}

#endif //BTER_SEQUENTIAL_BTERPHASES_H
