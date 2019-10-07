#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <dist/bicop_BB1_log.hpp>
#include <dist/bicop_survival_BB1_log.hpp>
#include <dist/bicop_r90_BB1_log.hpp>
#include <dist/bicop_r270_BB1_log.hpp>

#include <iostream>

TEST(Copula_density, BB1_copula_v) {
    using stan::math::var;

    double v_val = 0.8;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_BB1_log<false>(0.1, v, 0.5, 1.5);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,-1.679066);
        EXPECT_FLOAT_EQ(lp1adj,-4.55145);
}


TEST(Copula_density, BB1_copula_u) {
    using stan::math::var;

    double v_val = 0.05;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_BB1_log<false>(v, 0.7, 0.5, 1.5);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,-1.838721);
        EXPECT_FLOAT_EQ(lp1adj,16.320076);
}


TEST(Copula_density, BB1_copula_theta) {
    using stan::math::var;

    double theta_val = 0.5;
        var theta(theta_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_BB1_log<false>(0.1, 0.5, theta, 2.);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = theta.adj();


        EXPECT_FLOAT_EQ(lp1val,-1.048868);
        EXPECT_FLOAT_EQ(lp1adj,-1.683259);
}


TEST(Copula_density, BB1_copula_delta) {
    using stan::math::var;

    double theta_val = 2;
        var theta(theta_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_BB1_log<false>(0.1, 0.5, 0.5, theta);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = theta.adj();


        EXPECT_FLOAT_EQ(lp1val,-1.048868);
        EXPECT_FLOAT_EQ(lp1adj,-1.176625);
}

TEST(Copula_density, survival_BB1_copula_u) {
    using stan::math::var;

    double v_val = 0.05;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_survival_BB1_log<false>(v, 0.7, 0.5, 1.5);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,-1.546346);
        EXPECT_FLOAT_EQ(lp1adj,11.630591);
}

TEST(Copula_density, survival_BB1_copula_v) {
    using stan::math::var;

    double v_val = 0.8;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_survival_BB1_log<false>(0.1, v, 0.5, 1.5);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,-1.562345);
        EXPECT_FLOAT_EQ(lp1adj,-5.118207);
}



TEST(Copula_density, survival_BB1_copula_theta) {
    using stan::math::var;

    double theta_val = 0.5;
        var theta(theta_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_survival_BB1_log<false>(0.1, 0.5, theta, 2.);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = theta.adj();


        EXPECT_FLOAT_EQ(lp1val,-1.088868);
        EXPECT_FLOAT_EQ(lp1adj,-0.4009556);
}


TEST(Copula_density, survival_BB1_copula_delta) {
    using stan::math::var;

    double theta_val = 2;
        var theta(theta_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_survival_BB1_log<false>(0.1, 0.5, 0.5, theta);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = theta.adj();


        EXPECT_FLOAT_EQ(lp1val,-1.088868);
        EXPECT_FLOAT_EQ(lp1adj,-1.446326);
}

TEST(Copula_density, r90_BB1_copula_u) {
    using stan::math::var;

    double v_val = 0.05;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_r90_BB1_log<false>(v, 0.7, -0.5, -1.5);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,-0.06685771);
        EXPECT_FLOAT_EQ(lp1adj,9.9650469);
}

TEST(Copula_density, r90_BB1_copula_v) {
    using stan::math::var;

    double v_val = 0.8;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_r90_BB1_log<false>(0.1, v, -0.5, -1.5);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,0.6518384);
        EXPECT_FLOAT_EQ(lp1adj,4.1926026);
}



TEST(Copula_density, r90_BB1_copula_theta) {
    using stan::math::var;

    double theta_val = -0.5;
        var theta(theta_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_r90_BB1_log<false>(0.5, 0.5, theta, -2.);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = theta.adj();


        EXPECT_FLOAT_EQ(lp1val,0.5905375);
        EXPECT_FLOAT_EQ(lp1adj,-0.32509);
}


TEST(Copula_density, r90_BB1_copula_delta) {
    using stan::math::var;

    double theta_val = -2.;
        var theta(theta_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_r90_BB1_log<false>(0.1, 0.5, -0.5, theta);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = theta.adj();


        EXPECT_FLOAT_EQ(lp1val,-1.088868);
        EXPECT_FLOAT_EQ(lp1adj,1.44594);
}

TEST(Copula_density, r270_BB1_copula_u) {
    using stan::math::var;

    double v_val = 0.05;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_r270_BB1_log<false>(v, 0.7, -0.5, -1.5);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,-0.11505394);
        EXPECT_FLOAT_EQ(lp1adj,11.944067);
}

TEST(Copula_density, r270_BB1_copula_v) {
    using stan::math::var;

    double v_val = 0.8;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_r270_BB1_log<false>(0.1, v, -0.5, -1.5);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,0.6750373);
        EXPECT_FLOAT_EQ(lp1adj,4.3265877);
}



TEST(Copula_density, r270_BB1_copula_theta) {
    using stan::math::var;

    double theta_val = -0.5;
        var theta(theta_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_r270_BB1_log<false>(0.1, 0.5, theta, -2.);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = theta.adj();


        EXPECT_FLOAT_EQ(lp1val,-1.048868);
        EXPECT_FLOAT_EQ(lp1adj,1.682212);
}


TEST(Copula_density, r270_BB1_copula_delta) {
    using stan::math::var;

    double theta_val = -2;
        var theta(theta_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_r270_BB1_log<false>(0.1, 0.5, -0.5, theta);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = theta.adj();


        EXPECT_FLOAT_EQ(lp1val,-1.048868);
        EXPECT_FLOAT_EQ(lp1adj,1.176264);
}
//TEST(Copula_density, BB1_copula_v_50) {
//    using stan::math::var;
//
//    double v_val = 0.5;
//        var v(v_val);
//        var lp1(0.0);
//        lp1 += vifcopula::bicop_BB1_log<false>(0.1, v, 50);
//        double lp1val = lp1.val();
//
//        lp1.grad();
//        double lp1adj = v.adj();
//
//        EXPECT_FLOAT_EQ(lp1val,-75.84692);
//        EXPECT_FLOAT_EQ(lp1adj,-102);
//}
//
//TEST(Copula_density, BB1_copula_theta_50) {
//    using stan::math::var;
//
//    double theta_val = 50;
//        var theta(theta_val);
//        var lp1(0.0);
//        lp1 += vifcopula::bicop_BB1_log<false>(0.1, 0.5, theta);
//        double lp1val = lp1.val();
//
//        lp1.grad();
//        double lp1adj = theta.adj();
//
//
//        EXPECT_FLOAT_EQ(lp1adj,-1.58983);
//}
