#ifndef BTER_SEQUENTIAL_BTERPHASES_H
#define BTER_SEQUENTIAL_BTERPHASES_H

namespace BTERSeq {

    struct BterSamples {
        double* r;
        double s1;
        double s2;
    };

    class BterPhases {
        BterSamples computeSamples(double *wg, double *wd, int dmax);

        double randomUnified(int from, int to);

        double *randomSample(double *wg, double s1);
    };
}

#endif //BTER_SEQUENTIAL_BTERPHASES_H
