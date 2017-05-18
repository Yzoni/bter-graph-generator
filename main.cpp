#include <iostream>

#include "BterPhases.h"

#define DMAX 7

using namespace bter;

void setupLogger() {

}

extern void cudaTestFunction();

int main() {

    setupLogger();

    cudaTestFunction();

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

    BTERSetup bterSetup(nd, cd, &beta, dmax, &bterSetupResult);
    bterSetup.run();

    // Phases
    BterPhases bterPhases(&bterSetupResult, dmax, nd, cd);
    bterPhases.computeSamples();

    int *phase_one_i = new int[bterPhases.bterSamples.s1];
    int *phase_one_j = new int[bterPhases.bterSamples.s1];
    bterPhases.phaseTwo(phase_one_i, phase_one_j);

    int *phase_two_i = new int[bterPhases.bterSamples.s2];
    int *phase_two_j = new int[bterPhases.bterSamples.s2];
    bterPhases.phaseTwo(phase_two_i, phase_two_j);

//    std::chrono::time_point<std::chrono::system_clock> start, end;
//    start = std::chrono::system_clock::now();
//    std::cout << "f(42) = " << fibonacci(42) << '\n';
//    end = std::chrono::system_clock::now();
//
//    std::chrono::duration<double> elapsed_seconds = end-start;
//    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
//
//    std::cout << "finished computation at " << std::ctime(&end_time)
//              << "elapsed time: " << elapsed_seconds. << "s\n";
    return 0;
}

