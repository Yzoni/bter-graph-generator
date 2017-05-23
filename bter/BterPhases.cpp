#include <numeric>
#include <functional>
#include <cmath>
#include <random>
#include <iostream>
#include "BterPhases.h"


namespace bter {

    BterPhases::BterPhases(BTERSetupResult *bterSetupResult, int dmax, double *nd, double *cd) :
            bterSetupResult(bterSetupResult), dmax(dmax), nd(nd), cd(cd), bterSamples(0, 0) {

    }

    double BterPhases::randomUnified(int from, int to) {
        std::random_device rand_dev;
        std::mt19937 generator(rand_dev());
        std::uniform_real_distribution<double> distr(from, to);
        return distr(generator);
    }

    void BterPhases::randomSample(double *wg, int nsamples, int *binindices) {
        int last_element = dmax - 1; // hack

        double wg_sum = std::accumulate(wg, &wg[dmax], 0.0, std::plus<double>());
        for (int i = 1; i < dmax; ++i) {
            wg[i] = wg[i] * nsamples / wg_sum;
            if (wg[i] == 0) {
                last_element = i;
                break;
            }
        }

        // Generate bins
        std::vector<double> partial_sum_wg(last_element + 2);
        partial_sum_wg[0] = 0;
        std::partial_sum(wg, &wg[last_element], partial_sum_wg.begin() + 1);

        // Divide by total
        for (int c = 0; c < last_element + 1; ++c) {
            partial_sum_wg[c] = partial_sum_wg[c] / partial_sum_wg[last_element];
        }

        // Assign uniform random values to bins
        double r;
        for (int i = 0; i < nsamples; ++i) {
            r = randomUnified(0, 1);
            for (int j = 0; j < last_element + 1; ++j) {
                if (r > partial_sum_wg[j]) {
                    ++binindices[i];
                }
            }
        }
    }

    void BterPhases::computeSamples() {
        double wg_sum = std::accumulate(bterSetupResult->wg, &bterSetupResult->wg[dmax], 0.0, std::plus<double>());
        double wd_sum = std::accumulate(bterSetupResult->wd, &bterSetupResult->wd[dmax], 0.0, std::plus<double>());
        double w = wg_sum + wd_sum;

        int nsmp = (int) std::round(w);

        double *r = new double[nsmp];
        for (int i = 0; i < nsmp; ++i) {
            r[i] = randomUnified(0, 1);
        }

        double t = (wg_sum / w);

        for (int i = 0; i < nsmp; ++i) {
            if (r[i] < t) {
                ++bterSamples.s1;
            }
        }
        bterSamples.s2 = nsmp - bterSamples.s1;

        delete[] r;
    }

    void BterPhases::phaseOnePrepare(int *group_sample, double *block_b, double *block_i, double *block_n) {
        // Get group samples
        randomSample(bterSetupResult->wg, bterSamples.s1, group_sample);

        int k, sample;
        for (k = 0; k < bterSamples.s1; ++k) {
            sample = group_sample[k] - 1;
            block_b[k] = bterSetupResult->bg[sample]; // Index of first node in group g
            block_i[k] = bterSetupResult->ig[sample]; // Number of affinity block in group g
            block_n[k] = bterSetupResult->ng[sample]; // Number of of nodes in group g
        }
    }

    void BterPhases::phaseOneSeq(int *phase_one_i, int *phase_one_j) {

        int *group_sample = new int[bterSamples.s1]();
        double *block_b = new double[bterSamples.s1];
        double *block_i = new double[bterSamples.s1];
        double *block_n = new double[bterSamples.s1];

        phaseOnePrepare(group_sample, block_b, block_i, block_n);

        int k, shift;
        for (k = 0; k < bterSamples.s1; ++k) {
            // Compute block and compute its offset
            shift = (int) std::round(block_i[k] + floor(randomUnified(0, 1) * block_b[k]) * block_n[k]);

            // Choose first node
            phase_one_i[k] = (int) std::round(floor(randomUnified(0, 1) * block_n[k]) + shift);

            // Choose second node
            // "Without replacement"
            phase_one_j[k] = (int) std::round(floor(randomUnified(0, 1) * (block_n[k] - 1)) + shift);

            // Remove loops
            if (phase_one_j[k] >= phase_one_i[k]) {
                ++phase_one_j[k];
            }
        }

        delete[] group_sample;
        delete[] block_b;
        delete[] block_i;
        delete[] block_n;

    }

    void BterPhases::phaseOneGPU(int *phase_one_i, int *phase_one_j) {

        int *group_sample = new int[bterSamples.s1]();
        double *block_b = new double[bterSamples.s1];
        double *block_i = new double[bterSamples.s1];
        double *block_n = new double[bterSamples.s1];

        phaseOnePrepare(group_sample, block_b, block_i, block_n);

        cuda_wrapper_phase_one(phase_one_i, phase_one_j,
                               block_b, block_i, block_n,
                               bterSamples.s1);

        delete[] group_sample;
        delete[] block_b;
        delete[] block_i;
        delete[] block_n;
    }

    void BterPhases::phaseTwoSeq(int *phase_two_i, int *phase_two_j) {

        int i;
        double *id_bulk = new double[dmax];
        double *nd_bulk = new double[dmax];
        for (i = 0; i < dmax; ++i) {
            id_bulk[i] = bterSetupResult->id[i] + bterSetupResult->ndfill[i];
            nd_bulk[i] = nd[i] - bterSetupResult->ndfill[i];
        }

        phaseTwoNodeSeq(id_bulk, nd_bulk, phase_two_i);
        phaseTwoNodeSeq(id_bulk, nd_bulk, phase_two_j);

    }

    void BterPhases::phaseTwoGpu(int *phase_two_i, int *phase_one_j) {

        int i;
        double *id_bulk = new double[dmax];
        double *nd_bulk = new double[dmax];
        for (i = 0; i < dmax; ++i) {
            id_bulk[i] = bterSetupResult->id[i] + bterSetupResult->ndfill[i];
            nd_bulk[i] = nd[i] - bterSetupResult->ndfill[i];
        }

        phaseTwoNodeGpu(id_bulk, nd_bulk, phase_two_i);
        phaseTwoNodeGpu(id_bulk, nd_bulk, phase_two_j);
    }

    void BterPhases::phaseTwoNodePrepare(double *id_bulk, double *nd_bulk,
            double *phase_two_shift_fill, double *phase_two_sz_fill, double *phase_two_shift_bulk, double *phase_two_sz_bulk) {
        int *degree_sample = new int[bterSamples.s2]();

        // Excess degree sample
        randomSample(bterSetupResult->wd, bterSamples.s2, degree_sample);

        int i, sample;
        for (i = 0; i < bterSamples.s2; ++i) {
            sample = degree_sample[i] - 1;
            phase_two_shift_fill[i] = bterSetupResult->id[sample];
            phase_two_sz_fill[i] = bterSetupResult->ndfill[sample];
            phase_two_shift_bulk[i] = id_bulk[sample];
            pase_two_sz_bulk[i] = nd_bulk[sample];
        }
    }

    void BterPhases::phaseTwoNodeSeq(double *id_bulk, double *nd_bulk, int *phase_two) {

        double *phase_two_shift_fill = new double[bterSamples.s2];
        double *phase_two_sz_fill = new double[bterSamples.s2];
        double *phase_two_shift_bulk = new double[bterSamples.s2];
        double *pase_two_sz_bulk = new double[bterSamples.s2];

        phaseTwoNodePrepare(phase_two_shift_fill, phase_two_sz_fill, phase_two_shift_bulk, phase_two_sz_bulk);        

        // TODO GPU
        double phase_two_fill, phase_two_bulk, r;
        for (i = 0; i < bterSamples.s2; ++i) {
            phase_two_fill = phase_two_shift_fill[i] + floor(randomUnified(0, 1) * phase_two_sz_fill[i]);
            phase_two_bulk = phase_two_shift_bulk[i] + floor(randomUnified(0, 1) * pase_two_sz_bulk[i]);

            r = randomUnified(0, 1);
            if (r < phase_two_fill) {
                phase_two[i] = (int) std::round(phase_two_fill);
            } else {
                phase_two[i] = (int) std::round(phase_two_bulk);
            }
        }

        delete[] phase_two_shift_fill;
        delete[] phase_two_sz_fill;
        delete[] phase_two_shift_bulk;
        delete[] pase_two_sz_bulk;
    }


    void BterPhases::phaseTwoNodeGpu(double *id_bulk, double *nd_bulk, int *phase_two) {

        double *phase_two_shift_fill = new double[bterSamples.s2];
        double *phase_two_sz_fill = new double[bterSamples.s2];
        double *phase_two_shift_bulk = new double[bterSamples.s2];
        double *pase_two_sz_bulk = new double[bterSamples.s2];

        phaseTwoNodePrepare(phase_two_shift_fill, phase_two_sz_fill, phase_two_shift_bulk, phase_two_sz_bulk);        

        // TODO GPU
        double phase_two_fill, phase_two_bulk, r;
        for (i = 0; i < bterSamples.s2; ++i) {
            phase_two_fill = phase_two_shift_fill[i] + floor(randomUnified(0, 1) * phase_two_sz_fill[i]);
            phase_two_bulk = phase_two_shift_bulk[i] + floor(randomUnified(0, 1) * pase_two_sz_bulk[i]);

            r = randomUnified(0, 1);
            if (r < phase_two_fill) {
                phase_two[i] = (int) std::round(phase_two_fill);
            } else {
                phase_two[i] = (int) std::round(phase_two_bulk);
            }
        }

        delete[] phase_two_shift_fill;
        delete[] phase_two_sz_fill;
        delete[] phase_two_shift_bulk;
        delete[] pase_two_sz_bulk;
    }

    /*
     * Very naive implementation
     */
//    void BterPhases::removeLoopsPhaseTwo(int *i, int *j, int length, std::vector<int> *new_i, std::vector<int> *new_j) {
//        std::vector<int>::iterator it_new_i = new_i->begin();
//        std::vector<int>::iterator it_new_j = new_j->begin();
//
//        for (int k = 0; k < length; ++k) {
//            if (i[k] != j[k]) {
//                *it_new_i = i[k];
//                *it_new_j = i[k];
//
//                std::next(it_new_i);
//                std::next(it_new_j);
//            }
//        }
//
//
//        new_i->resize(new_i->size());
//        new_j->resize(new_j->size());
//    }
}
