#ifndef VIFCOPULA_DISTRIBUTION_BICOP_INDEPENDENCE_LOG_TEST_CPP
#define VIFCOPULA_DISTRIBUTION_BICOP_INDEPENDENCE_LOG_TEST_CPP

#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <dist/bicop_independence_log.cpp>


TEST(Copula_density, independence_copula) {
    using stan::math::var;

    double v_val = 0.1;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_independence_log<true>(0.1, v);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();
        std::cout << lp1val << " " << lp1adj << std::endl;

        EXPECT_FLOAT_EQ(lp1val,0.0);
        EXPECT_FLOAT_EQ(lp1adj,1.0);

}

#endif // VIFCOPULA_DISTRIBUTION_BICOP_INDEPENDENCE_LOG_TEST_CPP
