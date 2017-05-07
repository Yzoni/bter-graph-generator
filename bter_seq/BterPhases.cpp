#include <numeric>
#include <functional>
#include <cmath>
#include <random>
#include "BterPhases.h"


namespace BTERSeq {

    double BterPhases::randomUnified(int from, int to) {
        std::random_device rand_dev;
        std::mt19937 generator(rand_dev());
        std::uniform_int_distribution<double> distr(from, to);
        return distr(generator);
    }

    double *BterPhases::randomSample(double *wg, double s1) {

    }

    BterSamples BterPhases::computeSamples(double *wg, double *wd, int dmax) {
        double wg_sum = std::accumulate(wg, &wg[dmax - 1], 0, std::plus<int>());
        double wd_sum = std::accumulate(wd, &wd[dmax - 1], 0, std::plus<int>());

        double w = wg_sum + wd_sum;
        double nsmp = std::round(w);

        double *r = new double[nsmp];

        BterSamples bterSamples;
        bterSamples.r = r;

        for (int i = 0; i < nsmp; ++i) {
            r[i] = randomUnified(0, 1);
        }

        bterSamples.s1;
        double t = (wg_sum / w);
        for (int i = 0; i < nsmp; ++i) {
            if (r[i] < t) {
                ++bterSamples.s1;
            }
        }
        bterSamples.s2 = nsmp - bterSamples.s1;

        return bterSamples;
    }
}