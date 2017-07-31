#ifndef VIFCOPULA_HFUNC_JOE_HPP
#define VIFCOPULA_HFUNC_JOE_HPP

#include <stan/math.hpp>


namespace vifcopula
{

using namespace stan::math;
using namespace stan;

/**
 * The Conditional Distribution Function of a Bivariate Joe Copula
 * for the specified vector(s) u and the specified vector(s) v given parameter(s) theta.
 * u, v, or theta can each be either a scalar or a vector. Any vector inputs
 * must be the same length.
 *
 * <p>The result of h2 function is defined in VineCopula package.
 * @param u (Sequence of) scalar(s) in [0,1].
 * @param v (Sequence of) scalar(s) in [0,1].
 * @param theta (Sequence of) parameters for the Joe copula
 * @return The Conditional Distribution Function of a bivariate Joe Copula of u given v
 * @tparam T_u Underlying type of scalar in sequence.
 * @tparam T_v Underlying type of scalar in sequence.
 * @tparam T_theta Type of correlation parameter.
 */

template <typename T_u, typename T_v, typename T_theta>
inline
T_u hfunc_joe(const T_u& u, const T_v& v, const T_theta& theta)
{
    static const char* function("vifcopula::hfunc_joe");

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

    scalar_seq_view<T_u> u_vec(u);
    scalar_seq_view<T_v> v_vec(v);
    scalar_seq_view<T_theta> theta_vec(theta);
    size_t N = stan::max_size(u, v, theta);

    T_u u_cond(u);

    for (size_t n = 0; n < N; n++)
    {
        double u_val = value_of(u_vec[n]);
        double v_val = value_of(v_vec[n]);
        double theta_value = value_of(theta_vec[n]);


        double om_u = 1 - u_val;
        double om_v = 1 - v_val;
        double t_u = pow(om_u, theta_value);
        double t_v = pow(om_v, theta_value);
        double t_uv = t_u + t_v - t_u * t_v;

        u_cond[n] = pow(t_uv, 1/theta_value-1) * pow( om_v, theta_value -1 ) * (1 - t_u) ;

    }
    return u_cond;
}


}
#endif // VIFCOPULA_HFUNC_JOE_HPP


