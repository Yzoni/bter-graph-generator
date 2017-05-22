#include <valarray>

#ifndef BTER_SEQUENTIAL_BTERSEQ_H_H
#define BTER_SEQUENTIAL_BTERSEQ_H_H

namespace bter {

    struct BTERSetupResult {
        int *id, *ndprime;
        double *wd, *rdfill, *ndfill, *wg, *ig, *bg, *ng;
    };

    class BTERSetup {
    private:
        double *nd, *cd, *beta;
        int *id, *ndprime;
        double *wd, *rdfill, *ndfill, *wg, *ig, *bg, *ng;
        int dmax;

    public:
        int load_parameters();

        void compute_index_degree();

        void compute_degree_greater_then();

        void handle_degree_one();

        void run();

        BTERSetup(double *nd,
                  double *cd,
                  double *beta,
                  int dmax,
                  BTERSetupResult *bterSetupResult);
    };
}
#endif //BTER_SEQUENTIAL_BTERSEQ_H_H
