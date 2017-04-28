#ifndef VIFCOPULA_HFUNC_FRANK_HPP
#define VIFCOPULA_HFUNC_FRANK_HPP

#include <stan/math.hpp>


namespace vifcopula
{

using namespace stan::math;
using namespace stan;

/**
 * The Conditional Distribution Function of a Bivariate Frank Copula
 * for the specified vector(s) u and the specified vector(s) v given parameter(s) theta.
 * u, v, or theta can each be either a scalar or a vector. Any vector inputs
 * must be the same length.
 *
 * <p>The result of h2 function is defined in VineCopula package.
 * @param u (Sequence of) scalar(s) in [0,1].
 * @param v (Sequence of) scalar(s) in [0,1].
 * @param theta (Sequence of) parameters for the Frank copula
 * @return The Conditional Distribution Function of a bivariate Frank Copula of u given v
 * @tparam T_u Underlying type of scalar in sequence.
 * @tparam T_v Underlying type of scalar in sequence.
 * @tparam T_theta Type of correlation parameter.
 */

template <typename T_u, typename T_v, typename T_theta>
inline
T_u hfunc_frank(const T_u& u, const T_v& v, const T_theta& theta)
{
    static const char* function("vifcopula::hfunc_frank");

      using std::log;
      using std::pow;
      using std::exp;

//      check_bounded(function, "Random variable u", u,0,1);
//      check_bounded(function, "Random variable v", v,0,1);
//      check_bounded(function, "Correlation parameter", theta, -1,1);
    check_consistent_sizes(function,
                           "Random variable u", u,
                           "Random variable v", v,
                           "Scale parameter", theta);

    stan::VectorView<const T_u> u_vec(u);
    stan::VectorView<const T_v> v_vec(v);
    stan::VectorView<const T_theta> theta_vec(theta);
    size_t N = stan::max_size(u, v, theta);

    T_u u_cond(u);

    for (size_t n = 0; n < N; n++)
    {
        double u_val = value_of(u_vec[n]);
        double v_val = value_of(v_vec[n]);
        double theta_value = value_of(theta_vec[n]);
        double exp_theta = exp(theta_value);
        double exp_theta_u = pow( exp_theta, u_val);
        double exp_theta_v= pow( exp_theta, v_val);
        double exp_theta_uv= exp_theta_u * exp_theta_v;
        double exp_theta_up1 = exp_theta_u * exp_theta;
        double exp_theta_vp1 = exp_theta_v * exp_theta;

        u_cond[n] = - exp_theta * (exp_theta_v -1) / (exp_theta_uv - exp_theta_up1 - exp_theta_vp1 + exp_theta);

    }
    return u_cond;
}


}
#endif // VIFCOPULA_HFUNC_FRANK_HPP

