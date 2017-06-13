#include <numeric>
#include <functional>
#include <cmath>
#include <random>
#include <iostream>

#include "BterPhasesGpu.h"
#include "../spdlog/spdlog.h"

namespace spd = spdlog;

namespace bter {

    BterPhasesGpu::BterPhasesGpu(BTERSetupResult *bterSetupResult, int dmax, double *nd, double *cd) :
            bterSetupResult(bterSetupResult), dmax(dmax), nd(nd), cd(cd), bterSamples(0, 0) {

    }

    void BterPhasesGpu::randomSample(double *wg, int nsamples, int *binindices) {
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

        double *r = new double[nsamples];
        spd::get("logger")->info("Start randomSample() nsamples: {}", nsamples);
        cuda_wrapper_rand_array(nsamples, r);

        // Assign uniform random values to bins
        double ru;
        for (int i = 0; i < nsamples; ++i) {
            ru = r[i];
            for (int j = 0; j < last_element + 1; ++j) {
                if (ru > partial_sum_wg[j]) {
                    ++binindices[i];
                }
            }
        }

        delete[] r;
    }

    void BterPhasesGpu::computeSamples() {
        // Summed weight phase one
        double wg_sum = std::accumulate(bterSetupResult->wg, &bterSetupResult->wg[dmax], 0.0, std::plus<double>());

        // Summed weight phase two
        double wd_sum = std::accumulate(bterSetupResult->wd, &bterSetupResult->wd[dmax], 0.0, std::plus<double>());
        double w = wg_sum + wd_sum;

        int nsmp = (int) std::round(w);

        // Setup timer
        std::chrono::time_point<std::chrono::system_clock> start, end;
        std::chrono::duration<double> elapsed_seconds;

        spd::get("logger")->info("Start computeSamples() generate random, nsmp {}", nsmp);
        start = std::chrono::system_clock::now();
        double *r = new double[nsmp];
        cuda_wrapper_rand_array(nsmp, r);
        end = std::chrono::system_clock::now();
        elapsed_seconds = end - start;
        spd::get("logger")->info("Finished computeSamples() generate random, took {} seconds", elapsed_seconds.count());

        double t = (wg_sum / w);
        spd::get("logger")->info("Start computeSamples() t value {}", t);

        spd::get("logger")->info("Start computeSamples() samples s1");
        start = std::chrono::system_clock::now();
        for (int i = 0; i < nsmp; ++i) {
            if (r[i] < t) {
                ++bterSamples.s1;
            }
        }
        end = std::chrono::system_clock::now();
        elapsed_seconds = end - start;
        spd::get("logger")->info("computeSamples() samples s1, took {} seconds", elapsed_seconds.count());

        bterSamples.s2 = nsmp - bterSamples.s1;
        spd::get("logger")->info("computeSamples() samples nsmp: {}, s1: {}, s2: {}", nsmp, bterSamples.s1,
                                 bterSamples.s2);

        delete[] r;
    }

    void BterPhasesGpu::phaseOnePrepare(int *group_sample, double *block_b, double *block_i, double *block_n) {

        // Setup timer
        std::chrono::time_point<std::chrono::system_clock> start, end;
        std::chrono::duration<double> elapsed_seconds;

        // Get group samples
        spd::get("logger")->info("Start phaseOnePrepare() randomSample");
        start = std::chrono::system_clock::now();
        randomSample(bterSetupResult->wg, bterSamples.s1, group_sample);
        end = std::chrono::system_clock::now();
        elapsed_seconds = end - start;
        spd::get("logger")->info("Finished phaseOnePrepare() random sample, took {} seconds", elapsed_seconds.count());

        spd::get("logger")->info("Start phaseOnePrepare() loop");
        start = std::chrono::system_clock::now();
        int k, sample;
        for (k = 0; k < bterSamples.s1; ++k) {
            sample = group_sample[k] - 1;
            block_b[k] = bterSetupResult->bg[sample]; // Index of first node in group g
            block_i[k] = bterSetupResult->ig[sample]; // Number of affinity block in group g
            block_n[k] = bterSetupResult->ng[sample]; // Number of of nodes in group g
        }
        end = std::chrono::system_clock::now();
        elapsed_seconds = end - start;
        spd::get("logger")->info("Finished phaseOnePrepare() loop, took {} seconds", elapsed_seconds.count());

    }

    void BterPhasesGpu::phaseOne(int *phase_one_i, int *phase_one_j) {

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

    void BterPhasesGpu::phaseTwo(int *phase_two_i, int *phase_two_j) {

        double *id_bulk = new double[dmax];
        double *nd_bulk = new double[dmax];
        for (int i = 0; i < dmax; ++i) {
            id_bulk[i] = bterSetupResult->id[i] + bterSetupResult->ndfill[i];
            nd_bulk[i] = nd[i] - bterSetupResult->ndfill[i];
        }

        BterPhasesGpu::phaseTwoNode(id_bulk, nd_bulk, phase_two_i);
        BterPhasesGpu::phaseTwoNode(id_bulk, nd_bulk, phase_two_j);

        delete[] id_bulk;
        delete[] nd_bulk;
    }

    void BterPhasesGpu::phaseTwoNodePrepare(double *id_bulk, double *nd_bulk,
                                            double *phase_two_shift_fill, double *phase_two_sz_fill,
                                            double *phase_two_shift_bulk, double *phase_two_sz_bulk,
                                            double *phase_two_rd_fill) {
        int *degree_sample = new int[bterSamples.s2]();


        // Excess degree sample
        randomSample(bterSetupResult->wd, bterSamples.s2, degree_sample);

        int i, sample;
        for (i = 0; i < bterSamples.s2; ++i) {
            sample = degree_sample[i] - 1;
            phase_two_shift_fill[i] = bterSetupResult->id[sample];
            phase_two_sz_fill[i] = bterSetupResult->ndfill[sample];
            phase_two_shift_bulk[i] = id_bulk[sample];
            phase_two_sz_bulk[i] = nd_bulk[sample];

            phase_two_rd_fill[i] = bterSetupResult->rdfill[sample];
        }

        delete[] degree_sample;
    }

    void BterPhasesGpu::phaseTwoNode(double *id_bulk, double *nd_bulk, int *phase_two) {

        double *phase_two_shift_fill = new double[bterSamples.s2];
        double *phase_two_sz_fill = new double[bterSamples.s2];
        double *phase_two_shift_bulk = new double[bterSamples.s2];
        double *phase_two_sz_bulk = new double[bterSamples.s2];

        double *phase_two_rd_fill = new double[bterSamples.s2];

        // Get group samples

        phaseTwoNodePrepare(id_bulk, nd_bulk,
                            phase_two_shift_fill, phase_two_sz_fill,
                            phase_two_shift_bulk, phase_two_sz_bulk,
                            phase_two_rd_fill);

        cuda_wrapper_phase_two(phase_two,
                               phase_two_shift_fill, phase_two_sz_fill,
                               phase_two_shift_bulk, phase_two_sz_bulk,
                               phase_two_rd_fill,
                               bterSamples.s2);

        delete[] phase_two_shift_fill;
        delete[] phase_two_sz_fill;
        delete[] phase_two_shift_bulk;
        delete[] phase_two_sz_bulk;
    }

    /*
     * Very naive implementation
     */
//    void BterPhasesGpu::removeLoopsPhaseTwo(int *i, int *j, int length, std::vector<int> *new_i, std::vector<int> *new_j) {
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
