#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <dist/bicop_joe_log.cpp>

#include <iostream>

TEST(Copula_density, DISABLE_Joe_copula_v) {
    using stan::math::var;

    double v_val = 0.1;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_joe_log<false>(0.1, v, 2);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,0.5193628);
        EXPECT_FLOAT_EQ(lp1adj,-0.7530415);
}

TEST(Copula_density, DISABLE_Joe_copula_theta) {
    using stan::math::var;

    double theta_val = 10;
        var theta(theta_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_joe_log<false>(0.1, 0.5, theta);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = theta.adj();


        EXPECT_FLOAT_EQ(lp1val,-2.952879);
        EXPECT_FLOAT_EQ(lp1adj,-0.4829795);
}
