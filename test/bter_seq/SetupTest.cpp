#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "BterSetup.h"

#define DMAX 7

namespace BTERSeq {

    TEST(bter_seq_test, compute_index_degree) {

        double nd[]{2, 4, 7, 4, 1, 1, 1};

        double cd[]{0, 0.939066565437604, 0.327035150501861, 0.113901512792621, 0.038132502986394,
                    0.017869824183824, 0.005860223916729};

        double beta = 1;

        int dmax = DMAX;
        int id[DMAX] = {};
        double wd[DMAX] = {};
        double rdfill[DMAX] = {};
        double ndfill[DMAX] = {};
        double wg[DMAX] = {};
        double ig[DMAX] = {};
        double bg[DMAX] = {};
        double ng[DMAX]{};
        int ndprime[DMAX]{};


        BTERSetup bterSetup(nd, cd, &beta, id, wd, rdfill, ndfill, wg, ig, bg, ng, &dmax, ndprime);
        bterSetup.compute_index_degree();

        int correct[] = {19, 1, 5, 12, 16, 17, 18};
        ASSERT_THAT(id, testing::ContainerEq(correct));
    }

    TEST(bter_seq_test, compute_degree_greater_then) {

        double nd[]{2, 4, 7, 4, 1, 1, 1};

        double cd[]{0, 0.939066565437604, 0.327035150501861, 0.113901512792621, 0.038132502986394,
                    0.017869824183824, 0.005860223916729};

        double beta = 1;

        int dmax = DMAX;
        int id[DMAX] = {};
        double wd[DMAX] = {};
        double rdfill[DMAX] = {};
        double ndfill[DMAX] = {};
        double wg[DMAX] = {};
        double ig[DMAX] = {};
        double bg[DMAX] = {};
        double ng[DMAX]{};
        int ndprime[DMAX]{};


        BTERSetup bterSetup(nd, cd, &beta, id, wd, rdfill, ndfill, wg, ig, bg, ng, &dmax, ndprime);
        bterSetup.compute_index_degree();
        bterSetup.compute_degree_greater_then();

        int correct_ndprime[]{0, 14, 7, 3, 2, 1, 0};
        bterSetup.handle_degree_one();
        ASSERT_THAT(ndprime, testing::ContainerEq(correct_ndprime));
    }

    TEST(bter_seq_test, handle_degree_one) {

        double nd[]{2, 4, 7, 4, 1, 1, 1};

        double cd[]{0, 0.939066565437604, 0.327035150501861, 0.113901512792621, 0.038132502986394,
                    0.017869824183824, 0.005860223916729};

        double beta = 1;

        int dmax = DMAX;
        int id[DMAX] = {};
        double wd[DMAX] = {};
        double rdfill[DMAX] = {};
        double ndfill[DMAX] = {};
        double wg[DMAX] = {};
        double ig[DMAX] = {};
        double bg[DMAX] = {};
        double ng[DMAX]{};
        int ndprime[DMAX]{};


        BTERSetup bterSetup(nd, cd, &beta, id, wd, rdfill, ndfill, wg, ig, bg, ng, &dmax, ndprime);
        bterSetup.compute_index_degree();
        bterSetup.handle_degree_one();

        double correct_ndfill[]{2, 0, 0, 0, 0, 0, 0};
        double correct_wd[]{1, 0, 0, 0, 0, 0, 0};
        double correct_rdfill[]{1, 0, 0, 0, 0, 0, 0};
        ASSERT_THAT(ndfill, testing::ContainerEq(correct_ndfill));
        ASSERT_THAT(wd, testing::ContainerEq(correct_wd));
        ASSERT_THAT(rdfill, testing::ContainerEq(correct_rdfill));
    }


    TEST(bter_seq_test, run) {

        double nd[]{2, 4, 7, 4, 1, 1, 1};

        double cd[]{0, 0.939066565437604, 0.327035150501861, 0.113901512792621, 0.038132502986394,
                    0.017869824183824, 0.005860223916729};

        double beta = 1;

        int dmax = DMAX;
        int id[DMAX] = {};
        double wd[DMAX] = {};
        double rdfill[DMAX] = {};
        double ndfill[DMAX] = {};
        double wg[DMAX] = {};
        double ig[DMAX] = {};
        double bg[DMAX] = {};
        double ng[DMAX]{};
        int ndprime[DMAX]{};


        BTERSetup bterSetup(nd, cd, &beta, id, wd, rdfill, ndfill, wg, ig, bg, ng, &dmax, ndprime);
        bterSetup.run();

        // id,wd,ndfill,rdfill,ig,wg,bg,ng
        int corrent_id[]{19, 1, 5, 12, 16, 17, 18};
        double correct_wd[]{1,0.082952986928035,3.374227276603903,4.172538842141774,1.772888372257843,2.272888372257843,2.772888372257843};
        double correct_ndfill[]{2, 0, 2, 3, 1, 1, 1};
        double correct_rdfill[]{1, 0, 0.308656296120113, 0.694936723080458, 1, 1, 1};
        double correct_ig[]{1, 7, 15, 0, 0, 0, 0};
        double correct_wg[]{23.254653692022522, 14.014258268340225, 3.978514544364716, 0, 0, 0, 0};
        double correct_bg[]{2, 2, 1, 0, 0, 0, 0};
        double correct_ng[]{3, 4, 4, 0, 0, 0, 0};
        ASSERT_THAT(id, testing::ContainerEq(corrent_id));
        EXPECT_THAT(wd, testing::ContainerEq(correct_wd));
        EXPECT_THAT(ndfill, testing::ContainerEq(correct_ndfill));
        EXPECT_THAT(rdfill, testing::ContainerEq(correct_rdfill));
        EXPECT_THAT(ig, testing::ContainerEq(correct_ig));
        EXPECT_THAT(wg, testing::ContainerEq(correct_wg));
        EXPECT_THAT(bg, testing::ContainerEq(correct_bg));
        EXPECT_THAT(ng, testing::ContainerEq(correct_ng));
    }


}