#include <numeric>
#include <functional>
#include <cmath>
#include <random>
#include <iostream>
#include "BterPhases.h"


namespace BTERSeq {

    BterPhases::BterPhases(BTERSetupResult *bterSetupResult, int dmax) : bterSetupResult(bterSetupResult), dmax(dmax) {

    }

    double BterPhases::randomUnified(int from, int to) {
        std::random_device rand_dev;
        std::mt19937 generator(rand_dev());
        std::uniform_real_distribution<double> distr(from, to);
        return distr(generator);
    }

    void BterPhases::randomSample(double *wg, int nsamples, int *binindices) {
        double wg_sum = std::accumulate(wg, &wg[dmax - 1], 0, std::plus<int>());
        for (int i = 0; i < dmax; ++i) {
            wg[i] = wg[i] * nsamples / wg_sum;
        }

        // Generate bins
        std::vector<double> partial_sum_wg;
        partial_sum_wg[0] = 0;
        auto it = std::next(partial_sum_wg.begin(), 1);
        std::partial_sum(wg, &wg[dmax], it);

        // Assign uniform random values to bins
        double r;
        for (int i = 0; i < nsamples; ++i) {
            r = randomUnified(0, 1);
            for (int j = 0; i < dmax; ++i) {
                if (r < it[j]) {
                    ++binindices[j];
                    break;
                }
            }
        }
    }

    void BterPhases::computeNsmp(double *wg_sum, double *wd_sum, int *nsmp) {
        *wg_sum = std::accumulate(bterSetupResult->wg, &bterSetupResult->wg[dmax - 1], 0, std::plus<double>());
        *wd_sum = std::accumulate(bterSetupResult->wd, &bterSetupResult->wd[dmax - 1], 0, std::plus<double>());
        double w = *wg_sum + *wd_sum;
        *nsmp = (int) std::round(w);
    }

    void BterPhases::computeSamples(BterSamples *bterSamples) {
        double w = bterSamples->wg_sum + bterSamples->wd_sum;

        for (int i = 0; i < bterSamples->nsmp; ++i) {
            bterSamples->r[i] = randomUnified(0, 1);
        }

        bterSamples->s1;
        double t = (bterSamples->wg_sum / w);
        for (int i = 0; i < bterSamples->nsmp; ++i) {
            if (bterSamples->r[i] < t) {
                ++bterSamples->s1;
            }
        }
        bterSamples->s2 = bterSamples->nsmp - bterSamples->s1;
    }

    void BterPhases::phaseOne(int *phase_one) {
        int *group_sample = new int[bterSamples->s1];
        double *block_b = new double[bterSamples->s1];
        double *block_i = new double[bterSamples->s1];
        double *block_n = new double[bterSamples->s1];

        // Get group samples
        randomSample(bterSetupResult->wg, bterSamples->s1, group_sample);

        int k, sample;
        for (k = 0; k < bterSamples->s1; ++k) {
            sample = group_sample[k];
            block_b[k] = bterSetupResult->bg[sample]; // Index of first node in group g
            block_i[k] = bterSetupResult->ig[sample]; // Number of affinity block in group g
            block_n[k] = bterSetupResult->ng[sample]; // Number of of nodes in group g
        }

        // TODO GPU
        int i, j, shift;
        for (k = 0; k < bterSamples->s1; ++k) {
            // Compute block and compute its offset
            shift = (int) std::round(block_i[k] + floor(randomUnified(0, 1) * block_b[k]) * block_n[k]);

            // Choose first node
            i = (int) round(floor(randomUnified(0, 1) * block_n[k]) + shift);

            // Choose second node
            j = (int) round(floor(randomUnified(0, 1) * (block_n[k] - 1)) + shift);

            if (j >= i) {
                phase_one[k] = j;
            } else {
                phase_one[k] = j + 1;
            }
        }

        delete[] group_sample;
        delete[] block_b;
        delete[] block_i;
        delete[] block_n;
    }

    void BterPhases::phaseTwo(int *phase_two_i, int *phase_two_j) {

        int i;
        for (i = 0; i < dmax; ++i) {

        }

        int phase_two_fill;
        int phase_two_bulk;
        for (i = 0; i < bterSamples->s2; ++i) {

        }


    }

    void BterPhases::removeLoopsPhaseTwo() {

    }
}