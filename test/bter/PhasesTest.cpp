

#include <gtest/gtest.h>
#include <BterPhasesSeq.h>

#define DMAX 7

namespace bter {
    TEST(bter_seq_test, compute_samples) {

        double nd[]{5, 3, 5, 4, 2, 0, 1};

        double cd[]{0, 0.979080080307294, 0.362202768314448, 0.145664469924314, 0.038552033439531, 0.014016962863169,
                    0.005891345155566};

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
        BterPhasesSeq bterPhasesSeq(&bterSetupResult, dmax, nd, cd);

    bterPhasesSeq.

    computeSamples();

        std::cerr << "\ncompute samples: \n";
    std::cerr << "s1 Value: " << bterPhasesSeq.bterSamples.s1 <<
    std::endl;
    std::cerr << "s2 Value: " << bterPhasesSeq.bterSamples.s2 <<
    std::endl;
    }


    TEST(bter_seq_test, compute_phase_one) {

        double nd[]{5, 3, 5, 4, 2, 0, 1};

        double cd[]{0, 0.979080080307294, 0.362202768314448, 0.145664469924314, 0.038552033439531, 0.014016962863169,
                    0.005891345155566};

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
        BterPhasesSeq bterPhasesSeq(&bterSetupResult, dmax, nd, cd);

bterPhasesSeq.

computeSamples();

int *phase_one_i = new int[bterPhasesSeq.bterSamples.s1];
int *phase_one_j = new int[bterPhasesSeq.bterSamples.s1];
bterPhasesSeq.
phaseOne(phase_one_i, phase_one_j
);

        std::cerr << "\n BterPhases 1" << std::endl;
for (
int i = 0;
i<bterPhasesSeq.bterSamples.
s1;
++i)
            std::cerr << "[" << phase_one_i[i] << " - " << phase_one_j[i] << "] ";

        delete[] phase_one_i;
        delete[] phase_one_j;
    }

    TEST(bter_seq_test, compute_phase_two) {

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
        BterPhasesSeq bterPhasesSeq(&bterSetupResult, dmax, nd, cd);

bterPhasesSeq.

computeSamples();

int *phase_two_i = new int[bterPhasesSeq.bterSamples.s2];
int *phase_two_j = new int[bterPhasesSeq.bterSamples.s2];
bterPhasesSeq.
phaseTwo(phase_two_i, phase_two_j
);

        std::cerr << "\n BterPhases 2" << std::endl;
for (
int i = 0;
i<bterPhasesSeq.bterSamples.
s2;
++i)
            std::cerr << "[" << phase_two_i[i] << " - " << phase_two_j[i] << "] ";

        delete[] phase_two_i;
        delete[] phase_two_j;
    }
}
