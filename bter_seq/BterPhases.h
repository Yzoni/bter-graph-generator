#ifndef BTER_SEQUENTIAL_BTERPHASES_H
#define BTER_SEQUENTIAL_BTERPHASES_H

#include "BterSetup.h"

namespace BTERSeq {

    struct BterSamples {
        double *r;
        int s1;
        int s2;
        double wg_sum;
        double wd_sum;
        int nsmp;
    };

    class BterPhases {
    public:
        BterPhases(BTERSetupResult *bterSetupResult, int dmax);

        void computeSamples(BterSamples *bterSamples);

        void computeNsmp(double *wg_sum, double *wd_sum, int *nsmp);

    private:
        BterSamples *bterSamples;

        BTERSetupResult *bterSetupResult;

        int dmax;

        double randomUnified(int from, int to);

        void randomSample(double *wg, int nsamples, int *binindices);

        void phaseOne(int *phaseOne);

        void phaseTwo(int *phase_two_i, int *phase_two_j);

        void removeLoopsPhaseTwo();
    };
}

#endif //BTER_SEQUENTIAL_BTERPHASES_H
