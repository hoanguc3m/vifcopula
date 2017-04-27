#ifndef VIFCOPULA_HFUNC_INDEPENDENCE_HPP
#define VIFCOPULA_HFUNC_INDEPENDENCE_HPP

#include <stan/math.hpp>


namespace vifcopula
{

using namespace stan::math;
using namespace stan;

/**
 * The Conditional Distribution Function of a Bivariate Independence Copula
 * for the specified vector(s) u and the specified vector(s) v.
 * u, v can each be either a scalar or a vector. Any vector inputs
 * must be the same length.
 *
 * <p>The result of h function is defined in VineCopula package.
 * @param u (Sequence of) scalar(s) in [0,1].
 * @param v (Sequence of) scalar(s) in [0,1].
  * @return The Conditional Distribution Function of a bivariate independence Copula of u given v
 * @tparam T_u Underlying type of scalar in sequence.
 * @tparam T_v Underlying type of scalar in sequence.
*/

template <typename T_u, typename T_v>
inline
T_u hfunc_independence(const T_u& u, const T_v& v)
{
    static const char* function("vifcopula::hfunc_independence");

    if (!(stan::length(u)
            && stan::length(v)))
        return 0.0;

//      check_bounded(function, "Random variable u", u,0,1);
//      check_bounded(function, "Random variable v", v,0,1);
    check_consistent_sizes(function,
                           "Random variable u", u,
                           "Random variable v", v);

    T_u u_cond = u;

    return u_cond;
}



}
#endif // VIFCOPULA_HFUNC_INDEPENDENCE_HPP


