#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <transform/hfunc_independence.hpp>
#include <transform/hfunc_normal.hpp>
#include <transform/hfunc_student.hpp>

#include <iostream>

TEST(Hfunc, Hfunc_Indep) {

    double v = 0.1;
    double u = 0.1;

    double u_cond = vifcopula::hfunc_independence(u, v);
    EXPECT_FLOAT_EQ(u_cond,u);
}

TEST(Hfunc, Hfunc_GaussCop) {

    std::vector<double> v = {0.1};
    std::vector<double> u = {0.2};
    std::vector<double> u_cond_true = {0.1601363};
    double rho = 0.5;

    std::vector<double> u_cond = vifcopula::hfunc_normal(u, v, rho);
    //std::cout << u_cond[0] << std::endl;
    EXPECT_FLOAT_EQ(u_cond[0],u_cond_true[0]);
}

TEST(Hfunc, Hfunc_StuCop) {

    std::vector<double> v = {0.1};
    std::vector<double> u = {0.2};
    std::vector<double> u_cond_true = {0.139733};
    double rho = 0.5;
    double nu = 5;

    std::vector<double> u_cond = vifcopula::hfunc_student(u, v, rho, nu);
    //std::cout << u_cond[0] << std::endl;
    EXPECT_FLOAT_EQ(u_cond[0],u_cond_true[0]);
}
