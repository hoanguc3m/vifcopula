#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <dist/bicop_clayton_log.hpp>

#include <iostream>

TEST(Copula_density, Clayton_copula_v) {
    using stan::math::var;

    double v_val = 0.1;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_clayton_log<false>(0.1, v, 0.5);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,0.6239036);
        EXPECT_FLOAT_EQ(lp1adj,-3.121909);
}

TEST(Copula_density, Clayton_copula_theta) {
    using stan::math::var;

    double theta_val = 1;
        var theta(theta_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_clayton_log<false>(0.1, 0.5, theta);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = theta.adj();


        EXPECT_FLOAT_EQ(lp1val,-0.5090741);
        EXPECT_FLOAT_EQ(lp1adj,-0.7642303);
}

TEST(Copula_density, Clayton_copula_v_50) {
    using stan::math::var;

    double v_val = 0.5;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_clayton_log<false>(0.1, v, 50);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,-75.84692);
        EXPECT_FLOAT_EQ(lp1adj,-102);
}

TEST(Copula_density, Clayton_copula_theta_50) {
    using stan::math::var;

    double theta_val = 50;
        var theta(theta_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_clayton_log<false>(0.1, 0.5, theta);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = theta.adj();


        EXPECT_FLOAT_EQ(lp1adj,-1.58983);
}
//
//TEST(Copula_density, Clayton_copula_extreme_v) {
//    using stan::math::var;
//
//    double v_val = 1e-10;
//        var v(v_val);
//        var lp1(0.0);
//        lp1 += vifcopula::bicop_clayton_log<false>(1e-10, v, 50);
//        double lp1val = lp1.val();
//
//        lp1.grad();
//        double lp1adj = v.adj();
//
//        EXPECT_FLOAT_EQ(lp1val,25.05044);
//        EXPECT_FLOAT_EQ(lp1adj,0.03302823);
//}
//
//TEST(Copula_density, Clayton_copula_extreme_theta) {
//    using stan::math::var;
//
//    double theta_val = 50;
//        var theta(theta_val);
//        var lp1(0.0);
//        lp1 += vifcopula::bicop_clayton_log<false>(1e-10, 1e-10, theta);
//        double lp1val = lp1.val();
//
//        lp1.grad();
//        double lp1adj = theta.adj();
//
//
//        EXPECT_FLOAT_EQ(lp1val,25.05044);
//        EXPECT_FLOAT_EQ(lp1adj,0.03302823);
//}
