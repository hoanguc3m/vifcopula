#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <logBifcop.cpp>

#include <iostream>

TEST(Copula_density,  DISABLE_bi_copula_ind) {
    using stan::math::var;


    double v_val = 0.1;
    var v(v_val);
    var lp1(0.0);
    lp1 += vifcopula::logBifcop(0,0.1, v);
    double lp1val = lp1.val();

    lp1.grad();
    double lp1adj = v.adj();

    EXPECT_FLOAT_EQ(lp1val,0.0);
    EXPECT_FLOAT_EQ(lp1adj,1.0);

}

TEST(Copula_density,  DISABLE_bi_copula_gauss){
    using stan::math::var;

    double v_val = 0.1;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::logBifcop(1,0.1, v, 0.5);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        std::cout << lp1val << " " << lp1adj << std::endl;

        EXPECT_FLOAT_EQ(lp1val,0.6912992);
        EXPECT_FLOAT_EQ(lp1adj,-2.434119);
}


TEST(Copula_density, bi_copula_Student){
    using stan::math::var;

    double v_val = 0.1;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::logBifcop(2,0.1, v, 0.5, 4);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,0.8431354);
        EXPECT_FLOAT_EQ(lp1adj,-2.943159);
}

TEST(Copula_density, bi_copula_Clayton){
    using stan::math::var;

    double v_val = 0.1;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::logBifcop(3,0.1, v, 0.5);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,0.6239036);
        EXPECT_FLOAT_EQ(lp1adj,-3.121909);
}
