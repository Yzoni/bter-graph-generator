

#include <gtest/gtest.h>
#include <BterPhases.h>

#define DMAX 7

namespace BTERSeq {
    TEST(bter_seq_test, compute_samples) {

        double nd[]{2, 4, 7, 4, 1, 1, 1};

        double cd[]{0, 0.939066565437604, 0.327035150501861, 0.113901512792621, 0.038132502986394,
                    0.017869824183824, 0.005860223916729};

        double beta = 1;
        int dmax = DMAX;

        int id[DMAX]{};
        double wd[DMAX]{};
        double rdfill[DMAX]{};
        double ndfill[DMAX]{};
        double wg[DMAX]{};
        double ig[DMAX]{};
        double bg[DMAX]{};
        double ng[DMAX]{};
        int ndprime[DMAX]{};

        BTERSetupResult bterSetupResult{
                id, ndprime,
                wd, rdfill, ndfill, wg, ig, bg, ng
        };

        // Setup
        BTERSetup bterSetup(nd, cd, &beta, dmax, &bterSetupResult);
        bterSetup.run();

        // Phases
        BterPhases bterPhases(&bterSetupResult, dmax, nd, cd);

        bterPhases.computeSamples();

        std::cerr << "\ncompute samples: \n";
        std::cerr << "s1 Value: " << bterPhases.bterSamples.s1 << std::endl;
        std::cerr << "s2 Value: " << bterPhases.bterSamples.s2 << std::endl;
    }


    TEST(bter_seq_test, compute_phase_one) {

        double nd[]{2, 4, 7, 4, 1, 1, 1};

        double cd[]{0, 0.939066565437604, 0.327035150501861, 0.113901512792621, 0.038132502986394,
                    0.017869824183824, 0.005860223916729};

        double beta = 1;
        int dmax = DMAX;

        int id[DMAX]{};
        double wd[DMAX]{};
        double rdfill[DMAX]{};
        double ndfill[DMAX]{};
        double wg[DMAX]{};
        double ig[DMAX]{};
        double bg[DMAX]{};
        double ng[DMAX]{};
        int ndprime[DMAX]{};

        BTERSetupResult bterSetupResult{
                id, ndprime,
                wd, rdfill, ndfill, wg, ig, bg, ng
        };

        // Setup
        BTERSetup bterSetup(nd, cd, &beta, dmax, &bterSetupResult);
        bterSetup.run();

        // Phases
        BterPhases bterPhases(&bterSetupResult, dmax, nd, cd);

        bterPhases.computeSamples();

        int *phase_one_i = new int[bterPhases.bterSamples.s1];
        int *phase_one_j = new int[bterPhases.bterSamples.s2];
        bterPhases.phaseOne(phase_one_i, phase_one_j);

        delete[] phase_one_i;
        delete[] phase_one_j;
    }
}
