#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <dist/bicop_student_log.cpp>

#include <iostream>

TEST(Copula_density, DISABLE_Student_copula_v) {
    using stan::math::var;

    double v_val = 0.1;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_student_log<false>(0.1, v, 0.5, 4);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,0.8431354);
        EXPECT_FLOAT_EQ(lp1adj,-2.943159);
}

TEST(Copula_density, DISABLE_Student_copula_rho) {
    using stan::math::var;

    double rho_val = 0.5;
        var rho(rho_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_student_log<false>(0.1, 0.5, rho, 6);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = rho.adj();


        EXPECT_FLOAT_EQ(lp1val,-0.2500733);
        EXPECT_FLOAT_EQ(lp1adj,-1.015308);
}

TEST(Copula_density, DISABLE_Student_copula_nu) {
    using stan::math::var;

    double nu_val = 5;
        var nu(nu_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_student_log<false>(0.5, 0.5, 0.8, nu);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = nu.adj();


        std::cout << lp1val << " " << lp1adj << std::endl;

        EXPECT_FLOAT_EQ(lp1val,0.228148);
        EXPECT_FLOAT_EQ(lp1adj,0.005838773);
}

