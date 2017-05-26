#ifndef BTER_BTERPHASESSEQ_H
#define BTER_BTERPHASESSEQ_H

#include "BterSetup.h"

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

    class BterPhasesSeq {
    public:
        BterPhasesSeq(BTERSetupResult *bterSetupResult, int dmax, double *nd, double *cd);

        void computeSamples();

        void phaseOnePrepare(int *group_sample, double *block_b, double *block_i, double *block_n);

        void phaseOne(int *phase_one_i, int *phase_one_j);

        void phaseTwoNodePrepare(double *id_bulk, double *nd_bulk,
                                 double *phase_two_shift_fill, double *phase_two_sz_fill,
                                 double *phase_two_shift_bulk, double *phase_two_sz_bulk);

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

#endif //BTER_BTERPHASESSEQ_H
