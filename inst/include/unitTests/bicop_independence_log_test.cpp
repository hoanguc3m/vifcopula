#ifndef VIFCOPULA_DISTRIBUTION_BICOP_INDEPENDENCE_LOG_TEST_CPP
#define VIFCOPULA_DISTRIBUTION_BICOP_INDEPENDENCE_LOG_TEST_CPP

#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <dist/bicop_independence_log.cpp>


TEST(Copula_density, independence_copula) {
    using stan::math::var;
    for (double v_val = -.5; v_val = 0; v_val = 0.5) {
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_independence_log<true>(0, v);
        double lp1val = lp1.val();
        stan::math::grad(lp1.vi_);
        double lp1adj = lp1.adj();

        EXPECT_FLOAT_EQ(lp1val,0.0);
        EXPECT_FLOAT_EQ(lp1adj,1);
    }
}

#endif // VIFCOPULA_DISTRIBUTION_BICOP_INDEPENDENCE_LOG_TEST_CPP
