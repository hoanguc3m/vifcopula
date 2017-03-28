#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <dist/bicop_gumbel_log.hpp>

#include <iostream>

TEST(Copula_density, Gumbel_copula_v) {
    using stan::math::var;

    double v_val = 0.3;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_gumbel_log(0.1, v, 2);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,0.3437033);
        EXPECT_FLOAT_EQ(lp1adj,-3.203365);
}

TEST(Copula_density, Gumbel_copula_theta) {
    using stan::math::var;

    double theta_val = 10;
        var theta(theta_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_gumbel_log<false>(0.1, 0.5, theta);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = theta.adj();


        EXPECT_FLOAT_EQ(lp1val,-8.520775);
        EXPECT_FLOAT_EQ(lp1adj,-1.112054);
}

TEST(Copula_density, Gumbel_copula_v_50) {
    using stan::math::var;

    double v_val = 0.5;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_gumbel_log<false>(0.1, v, 50);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,-55.02987);
        EXPECT_FLOAT_EQ(lp1adj,-143.3841);
}

TEST(Copula_density, Gumbel_copula_theta_50) {
    using stan::math::var;

    double theta_val = 50;
        var theta(theta_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_gumbel_log<false>(0.1, 0.5, theta);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = theta.adj();


        EXPECT_FLOAT_EQ(lp1adj,-1.181053);
}
