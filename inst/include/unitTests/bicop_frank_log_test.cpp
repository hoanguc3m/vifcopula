#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <dist/bicop_frank_log.cpp>

#include <iostream>

TEST(Copula_density, Frank_copula_v) {
    using stan::math::var;

    double v_val = 0.1;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_frank_log<false>(0.1, v, 0.5);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        std::cout << lp1val << " " << lp1adj << std::endl;

        EXPECT_FLOAT_EQ(lp1val,0.1517319);
        EXPECT_FLOAT_EQ(lp1adj,-0.3813779);
}

TEST(Copula_density, Frank_copula_theta) {
    using stan::math::var;

    double theta_val = 1;
        var theta(theta_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_frank_log<false>(0.1, 0.5, theta);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = theta.adj();


        std::cout << lp1val << " " << lp1adj << std::endl;

        EXPECT_FLOAT_EQ(lp1val,-0.01920139);
        EXPECT_FLOAT_EQ(lp1adj,-0.0384675);
}
