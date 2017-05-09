

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

        BterSamples bterSamples;
        bterPhases.computeSamples(&bterSamples);

        std::cout << "compute samples: \n" << std::endl;
        std::cout << "s1 Value: " << bterSamples.s1 << std::endl;
        std::cout << "s2 Value: " << bterSamples.s2 << std::endl;
    }

}
