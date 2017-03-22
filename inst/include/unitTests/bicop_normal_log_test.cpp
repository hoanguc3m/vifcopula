#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <dist/bicop_normal_log.hpp>

#include <iostream>

TEST(Copula_density, DISABLE_Normal_copula_v) {
    using stan::math::var;

    double v_val = 0.1;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_normal_log<false>(0.1, v, 0.5);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,0.6912992);
        EXPECT_FLOAT_EQ(lp1adj,-2.434119);
}

TEST(Copula_density, DISABLE_Normal_copula_rho) {
    using stan::math::var;

    double rho_val = 0.5;
        var rho(rho_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_normal_log<false>(0.1, 0.5, rho);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = rho.adj();


        EXPECT_FLOAT_EQ(lp1val,-0.129888);
        EXPECT_FLOAT_EQ(lp1adj,-0.7932217);
}

