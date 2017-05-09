

#define DMAX 7

namespace BTERSeq {
//    TEST(bter_seq_test, compute_samples) {
//
//        double nd[]{2, 4, 7, 4, 1, 1, 1};
//
//        double cd[]{0, 0.939066565437604, 0.327035150501861, 0.113901512792621, 0.038132502986394,
//                    0.017869824183824, 0.005860223916729};
//
//        double beta = 1;
//        int dmax = DMAX;
//
//        int id[DMAX]{};
//        double wd[DMAX]{};
//        double rdfill[DMAX]{};
//        double ndfill[DMAX]{};
//        double wg[DMAX]{};
//        double ig[DMAX]{};
//        double bg[DMAX]{};
//        double ng[DMAX]{};
//        int ndprime[DMAX]{};
//
//        BTERSetupResult bterSetupResult{
//                id, ndprime,
//                wd, rdfill, ndfill, wg, ig, bg, ng
//        };
//
//        BTERSetup bterSetup(nd, cd, &beta, dmax, &bterSetupResult);
//        bterSetup.run();
//
//        BterPhases bterPhases(bterSetupResult, dmax);
//
//        double wg_sum, wd_sum;
//        int nsmp;
//
//        bterPhases.computeNsmp(wg, wd, &wg_sum, &wd_sum, &nsmp);
//
//        BterSamples *bterSamples;
//        bterSamples->wg_sum = wg_sum;
//        bterSamples->wd_sum = wd_sum;
//        bterSamples->nsmp = nsmp;
//        bterSamples->r = new double[nsmp];
//
//        bterPhases.computeSamples(bterSamples);
//
//        std::cout << "compute samples: \n" << std::endl;
//        std::cout << "s1 Value: " << bterSamples->s1 << std::endl;
//        std::cout << "s2 Value: " << bterSamples->s2 << std::endl;
//        for (int i = 0; i < dmax - 1; ++i) {
//            std::cout << bterSamples->r[i] << " ";
//        }
//
//        delete[] bterSamples->r;
//    }
//
//    TEST(bter_seq_test, random_samples) {
//
//        double nd[]{2, 4, 7, 4, 1, 1, 1};
//
//        double cd[]{0, 0.939066565437604, 0.327035150501861, 0.113901512792621, 0.038132502986394,
//                    0.017869824183824, 0.005860223916729};
//
//        double beta = 1;
//
//        int dmax = DMAX;
//        int id[DMAX] = {};
//        double wd[DMAX] = {};
//        double rdfill[DMAX] = {};
//        double ndfill[DMAX] = {};
//        double wg[DMAX] = {};
//        double ig[DMAX] = {};
//        double bg[DMAX] = {};
//        double ng[DMAX]{};
//        int ndprime[DMAX]{};
//
//        BTERSetup bterSetup(nd, cd, &beta, id, wd, rdfill, ndfill, wg, ig, bg, ng, &dmax, ndprime);
//        bterSetup.run();
//
//        BterPhases bterPhases;
//        BterSamples bterSamples;
//        bterPhases.computeSamples(&bterSamples);
//
//        int *samples = new int[bterSamples.s1];
//        bterPhases.randomSample(wg, bterSamples.s1, dmax, samples);
//
//        delete[] bterSamples.r;
//        delete[] samples;
//    }
}
