#include <numeric>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <spdlog/spdlog.h>

#include "BterSetup.h"

namespace bter {

    BTERSetup::BTERSetup(
            double *nd,
            double *cd,
            double *beta,
            int dmax,
            BTERSetupResult *bterSetupResult
    ) : nd(nd), cd(cd),
        beta(beta),
        dmax(dmax),
        id(bterSetupResult->id),
        wd(bterSetupResult->wd),
        rdfill(bterSetupResult->rdfill),
        ndfill(bterSetupResult->ndfill),
        wg(bterSetupResult->wg),
        ig(bterSetupResult->ig),
        bg(bterSetupResult->bg),
        ng(bterSetupResult->ng),
        ndprime(bterSetupResult->ndprime) {
    }

    /**
     * Everything but the top 2 elements, after the inverse cumsum,
%    * with a 0 on the first and last element appended.
     * @param ndprime
     */
    void BTERSetup::compute_degree_greater_then() {
        std::vector<int> temp1(nd, &nd[dmax]);
        std::vector<int> temp2(nd, &nd[dmax]);

        std::partial_sum(temp1.rbegin(), rend(temp1), rbegin(temp2));
        auto it = std::next(temp2.begin(), 2);
        ndprime[0] = 0;
        std::copy(it, end(temp2), ndprime + 1);
        ndprime[dmax - 1] = 0;
    }

    /**
     * Index of first node for each degree.
     * Degree 1 vertices are numbered last.
     */
    void BTERSetup::compute_index_degree() {
        std::partial_sum(&nd[1], &nd[dmax], &id[0]);
        std::rotate(&id[0], &id[dmax - 2], &id[dmax]);
        for (int i = 0; i < dmax; ++i) ++id[i];
        id[1] = 1;
    }

    /**
     * Multiply with blowup factor
     */
    void BTERSetup::handle_degree_one() {
        ndfill[0] = nd[0] * *beta;
        rdfill[0] = nd[0] * 0.5;
        wd[0] = 0.5 * nd[0];
    }

    void BTERSetup::run() {
        double nfillblk = 0;
        double wdfilltmp = 0;
        double intdeg = 0;
        double ndbulktmp = 0;
        int g = -1;
        double rho = 0;
        double wdbulktmp = 0;

        compute_index_degree();
        compute_degree_greater_then();
        handle_degree_one();
        for (int d = 1; d < dmax; ++d) {
            // Try to fill incomplete block from current group
            if (nfillblk > 0) {
                ndfill[d] = std::min(nfillblk, nd[d]);
                nfillblk = nfillblk - ndfill[d];
                wdfilltmp = 0.5 * ndfill[d] * ((d + 1) - intdeg);
            } else {
                ndfill[d] = 0;
                wdfilltmp = 0;
            };

            ndbulktmp = nd[d] - ndfill[d];

            // Create a new group for degree-d bulk nodes
            if (ndbulktmp > 0) {
                g = g + 1;
                ig[g] = id[d] + ndfill[d];
                bg[g] = std::ceil(ndbulktmp / (d + 2));
                ng[g] = d + 2;

                // Special handling for last group
                if ((bg[g] * (d + 2)) > (ndprime[d] + ndbulktmp)) {
                    if (bg[g] != 1) {
                        std::cerr << "Last group has more than 1 block" << std::endl;
                    }
                    ng[g] = ndprime[d] + ndbulktmp;
                }
                rho = std::pow(cd[d], (1.0 / 3.0));
                intdeg = (ng[g] - 1.0) * rho;
                wdbulktmp = 0.5 * ndbulktmp * (d + 1 - intdeg);

                // Calculate weight for picking the group in phase 1
                wg[g] = bg[g] * 0.5 * ng[g] * (ng[g] - 1) * std::log(1.0 / (1 - rho));
                nfillblk = bg[g] * ng[g] - ndbulktmp;
            } else {
                wdbulktmp = 0;
            }

            // Calculate excess degree for picking degree d in phase 2
            wd[d] = wdbulktmp + wdfilltmp;
            if (wd[d] > 0) {
                rdfill[d] = wdfilltmp / wd[d];
            } else {
                rdfill[d] = 0;
            }
        }
    }

    int BTERSetup::load_parameters() {
        return 0;
    }
}