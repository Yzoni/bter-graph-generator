#ifndef BTER_SEQUENTIAL_BTERPHASES_H
#define BTER_SEQUENTIAL_BTERPHASES_H

#include "BterSetup.h"

void cuda_wrapper_phase_one(int *phase_one_i, int *phase_one_j,
                            double *block_b, double *block_i, double *block_n,
                            int length);

void cuda_wrapper_phase_two();

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

        void computeSamples();

        void phaseOnePrepare(int *group_sample, double *block_b, double *block_i, double *block_n);

        void phaseOneSeq(int *phase_one_i, int *phase_one_j);

        void phaseOneGPU(int *phase_one_i, int *phase_one_j);

        void phaseOne(int *phase_one_i, int *phase_one_j);

        void phaseTwo(int *phase_two_i, int *phase_two_j);

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