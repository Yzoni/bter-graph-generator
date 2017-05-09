#ifndef BTER_SEQUENTIAL_BTERPHASES_H
#define BTER_SEQUENTIAL_BTERPHASES_H

#include "BterSetup.h"

namespace BTERSeq {

    struct BterSamples {
        int s1, s2;
    };

    class BterPhases {
    public:
        BterPhases(BTERSetupResult *bterSetupResult, int dmax, double *nd, double *cd);

        void computeSamples(BterSamples *bterSamples);

    private:
        BterSamples *bterSamples;

        BTERSetupResult *bterSetupResult;

        int dmax;

        double *nd, *cd;

        double randomUnified(int from, int to);

        void randomSample(double *wg, int nsamples, int *binindices);

        void phaseOne(int *phaseOne);

        void phaseTwo(int *phase_two_i, int *phase_two_j);

        void removeLoopsPhaseTwo();

        void phaseTwoOneEdge(double *id_bulk, double *nd_bulk, int *phase_two);

        void phaseTwoNode(double *id_bulk, double *nd_bulk, int *phase_two);
    };
}

#endif //BTER_SEQUENTIAL_BTERPHASES_H
