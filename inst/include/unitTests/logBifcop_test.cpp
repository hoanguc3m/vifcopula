#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <logBifcop.cpp>

#include <iostream>

TEST(Copula_density, bi_copula_ind) {
    using stan::math::var;


        double v_val = 0.1;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::logBifcop(2,0.1, v,0.0,0.0);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();
        std::cout << lp1val << " " << lp1adj << std::endl;

        EXPECT_FLOAT_EQ(lp1val,0.0);
        EXPECT_FLOAT_EQ(lp1adj,1.0);

}

TEST(Copula_density, bi_copula_gauss){

}


