#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <dist/bicop_frank_log.hpp>
#include <dist/bicop_r90_frank_log.hpp>
#include <boost/random/additive_combine.hpp> // L'Ecuyer RNG
#include <iostream>


typedef boost::ecuyer1988 rng_t;

using stan::math::uniform_rng;
using std::vector;


TEST(Copula_density, Frank_copula_v) {
    using stan::math::var;

    double v_val = 0.1;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_frank_log<false>(0.1, v, 0.5);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

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


        EXPECT_FLOAT_EQ(lp1val,-0.01920139);
        EXPECT_FLOAT_EQ(lp1adj,-0.0384675);
}

TEST(Copula_density, Frank_copula_v_50) {
    using stan::math::var;

    double v_val = 0.5;
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_frank_log<false>(0.1, v, 50);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = v.adj();

        EXPECT_FLOAT_EQ(lp1val,-16.08798);
        EXPECT_FLOAT_EQ(lp1adj,-50);
}


TEST(Copula_density, Frank_copula_theta_50) {
    using stan::math::var;

    double theta_val = 100;
    var theta(theta_val);
    var lp1(0.0);

    // Set seed
    int seed  = 10;
    rng_t base_rng(seed);

    int t_max  = 11;

    vector<double> u(t_max);
    for (int i = 0; i < t_max; i++ )
    {
        u[i] = uniform_rng(0.0,1.0,base_rng);
        //std::cout << "u[i] " << u[i]  << std::endl;
    }



    vector<double> v(t_max);
    for (int i = 0; i < t_max; i++ )
    {
        //v[i] = (u[i] + uniform_rng(0.0,1.0,base_rng))/2;
        v[i] = uniform_rng(0.0,1.0,base_rng);
        //std::cout << "v[i] " << v[i]  << std::endl;
    }
    int copula_type = 5;

    lp1 += vifcopula::bicop_frank_log<false>(u, v, theta);
    double lp1val = lp1.val();

    lp1.grad();
    double lp1adj = theta.adj();
//    EXPECT_FLOAT_EQ(lp1val,-39.15095);
//    EXPECT_FLOAT_EQ(lp1adj,-0.8995);
}

TEST(Copula_density, Neg_Frank_copula_theta) {
    using stan::math::var;

    double theta_val = -1;
        var theta(theta_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_r90_frank_log<false>(0.1, 0.4, theta);
        double lp1val = lp1.val();

        lp1.grad();
        double lp1adj = theta.adj();


        EXPECT_FLOAT_EQ(lp1val,-0.1006427);
        EXPECT_FLOAT_EQ(lp1adj,0.12184646);
}
