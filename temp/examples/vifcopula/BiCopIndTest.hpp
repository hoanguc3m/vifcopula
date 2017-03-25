#ifndef VIFCOPULA_BiCopIndTest_HPP
#define VIFCOPULA_BiCopIndTest_HPP


/*
    Based on the BiCopIndTest function in VineCopula
    Author: Jeffrey Dissmann
    Object: test based on Kendall's of a bivariate independence copula
    @param u,v vectors of equal length with values in [0,1].
    @return P-value of the independence test.
    @references Genest, C. and A. C. Favre (2007). Everything you always wanted
    to know about copula modeling but were afraid to ask. Journal of Hydrologic
    Engineering, 12 (4), 347-368.
*/

#include <fastkendall.hpp>
#include <stan/math.hpp>


template <typename T_u, typename T_v>
double BiCopIndTest(const T_u& u, const T_v& v){

    static const char* function("vifcopula::BiCopIndTest");

    size_t t_max = stan::length(u);
    //size_t t_max2 = length(v);
    //stan::math::equal(function, "Number of vector elements",t_max, t_max2);

    double tau = fast_kendall(u, v, t_max);

    double f = sqrt((9 * t_max * (t_max - 1))/(2 * (2 * t_max + 5))) * abs(tau);
    boost::math::normal s;

    double p_value = 2 * (1 - pdf(s,f));
    return(p_value);

    }






#endif // VIFCOPULA_BiCopIndTest_HPP
