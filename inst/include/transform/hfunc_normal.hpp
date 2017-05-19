#ifndef VIFCOPULA_HFUNC_NORMAL_HPP
#define VIFCOPULA_HFUNC_NORMAL_HPP

#include <stan/math.hpp>


namespace vifcopula
{

using namespace stan::math;
using namespace stan;

/**
 * The Conditional Distribution Function of a Bivariate Copula
 * for the specified vector(s) u and the specified vector(s) v given correlation(s) rho.
 * u, v, or rho can each be either a scalar or a vector. Any vector inputs
 * must be the same length.
 *
 * <p>The result of h2 function is defined in VineCopula package.
 * @param u (Sequence of) scalar(s) in [0,1].
 * @param v (Sequence of) scalar(s) in [0,1].
 * @param rho (Sequence of) correlation parameters for the normal copula
 * @return The Conditional Distribution Function of a bivariate Gaussian Copula of u given v
 * @tparam T_u Underlying type of scalar in sequence.
 * @tparam T_v Underlying type of scalar in sequence.
 * @tparam T_rho Type of correlation parameter.
 */

template <typename T_u, typename T_v, typename T_rho>
inline
T_u hfunc_normal(const T_u& u, const T_v& v, const T_rho& rho)
{
    static const char* function("vifcopula::hfunc_normal");

    using stan::math::square;
    using stan::math::sqrt;
    using boost::math::normal;

//      check_bounded(function, "Random variable u", u,0,1);
//      check_bounded(function, "Random variable v", v,0,1);
//      check_bounded(function, "Correlation parameter", rho, -1,1);
    check_consistent_sizes(function,
                           "Random variable u", u,
                           "Random variable v", v,
                           "Scale parameter", rho);

    stan::VectorView<const T_u> u_vec(u);
    stan::VectorView<const T_v> v_vec(v);
    stan::VectorView<const T_rho> rho_vec(rho);
    size_t N = stan::max_size(u, v, rho);

    T_u u_cond(u);
    normal s;

    for (size_t n = 0; n < N; n++)
    {
        double inv_u_dbl = quantile(s,value_of(u_vec[n]));
        double inv_v_dbl = quantile(s,value_of(v_vec[n]));
        double rho_value = value_of(rho_vec[n]);
        double div_sigma = sqrt(1 - square(rho_value));

        u_cond[n] = cdf(s, (inv_u_dbl - rho_value * inv_v_dbl )/div_sigma );

    }
    return u_cond;
}


}
#endif // VIFCOPULA_HFUNC_NORMAL_HPP


