#include <valarray>

#ifndef BTER_SEQUENTIAL_BTERSEQ_H_H
#define BTER_SEQUENTIAL_BTERSEQ_H_H

namespace BTERSeq {

    class BTERSetup {
    private:
        double *nd;
        double *cd;
        double *beta;
        int *id;
        double *wd;
        double *rdfill;
        double *ndfill;
        double *wg;
        double *ig;
        double *bg;
        double *ng;
        int *dmax;
        int *ndprime;


    public:
        void compute_index_degree();

        void compute_degree_greater_then();

        void handle_degree_one();

        void run();

        BTERSetup(double *nd,
                  double *cd,
                  double *beta,
                  int *id,
                  double *wd,
                  double *rdfill,
                  double *ndfill,
                  double *wg,
                  double *ig,
                  double *bg,
                  double *ng,
                  int *dmax,
                  int *ndprime);
    };
}
#endif //BTER_SEQUENTIAL_BTERSEQ_H_H
