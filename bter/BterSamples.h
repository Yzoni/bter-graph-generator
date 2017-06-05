//
// Created by Yorick de Boer on 31/05/2017.
//

#ifndef BTER_PROJECT_BTERSAMPLES_H
#define BTER_PROJECT_BTERSAMPLES_H
namespace bter {

    struct BterSamples {
        int s1, s2;

        BterSamples(int s1, int s2) : s1(s1), s2(s2) {}
    };

    struct BterEdges {
        int *phase_one_i;
        int *phase_one_j;
        int *phase_two_i;
        int *phase_two_j;

        int bter_phase_one_size;
        int bter_phase_two_size;
    };
}
#endif //BTER_PROJECT_BTERSAMPLES_H
