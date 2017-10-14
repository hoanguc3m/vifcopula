#ifndef VIFCOPULA_DISTRIBUTION_BICOP_ROTATE_FRANK_LOG_HPP
#define VIFCOPULA_DISTRIBUTION_BICOP_ROTATE_FRANK_LOG_HPP

#include <stan/math.hpp>

namespace vifcopula {

using namespace stan::math;
using namespace stan;

    /**
     * The log of the bivariate Frank copula density (negative side) for the specified vector(s) u
     * and the specified vector(s) v given degree(s) of dependence theta.
     *  u, v, or theta can each be either a scalar or a vector. Any vector inputs
     * must be the same length.
     *
     * <p>The result log probability is defined to be the sum of the
     * log probabilities for each observation/observation/dependence.
     * @param u (Sequence of) scalar(s) in [0,1].
     * @param v (Sequence of) scalar(s) in [0,1].
     * @param theta (Sequence of) dependence parameters for the Frank copula
     * @return The log of the product of the densities.
     * @tparam T_u Underlying type of scalar in sequence.
     * @tparam T_v Underlying type of scalar in sequence.
     * @tparam T_theta Type of dependence parameter.
     */

    template <bool propto, typename T_u, typename T_v, typename T_theta>
    typename stan::return_type<T_u, T_v, T_theta>::type
    bicop_r90_frank_log(const T_u& u, const T_v& v, const T_theta& theta) {
      static const char* function("vifcopula::bicop_r90_frank_log");
      typedef typename stan::partials_return_type<T_u, T_v, T_theta>::type
        T_partials_return;

      using std::log;
      using std::pow;
      using std::exp;

      if (!(stan::length(u)
            && stan::length(v)
            && stan::length(theta)))
        return 0.0;

      T_partials_return logp(0.0);

//      check_bounded(function, "Random variable u", u,0,1);
//      check_bounded(function, "Random variable v", v,0,1);
//      check_positive(function, "dependence parameter", square(theta) ); // theta^2 > 0
      check_consistent_sizes(function,
                             "Random variable u", u,
                             "Random variable v", v,
                             "Scale parameter", theta);
      if (!include_summand<propto, T_u, T_v, T_theta>::value)
        return 0.0;

      operands_and_partials<T_u, T_v, T_theta> ops_partials(u, v, theta);

      scalar_seq_view<T_u> u_vec(u);
      scalar_seq_view<T_v> v_vec(v);
      scalar_seq_view<T_theta> theta_vec(theta);

      size_t N = stan::max_size(u, v, theta);

      stan::VectorBuilder<true, T_partials_return, T_theta> theta_value(length(theta));
      stan::VectorBuilder<true, T_partials_return, T_theta> exp_theta(length(theta));


      for (size_t i = 0; i < length(theta); i++) {
        theta_value[i] = value_of(theta_vec[i]);
        exp_theta[i] = exp(theta_value[i]);
      }


      for (size_t n = 0; n < N; n++) {
        const T_partials_return u_dbl = value_of(u_vec[n]);
        const T_partials_return v_dbl = value_of(v_vec[n]);
        const T_partials_return exp_theta_u = pow( exp_theta[n], u_dbl);
        const T_partials_return exp_theta_v= pow( exp_theta[n], v_dbl);
        const T_partials_return exp_theta_uv= exp_theta_u * exp_theta_v;
        const T_partials_return exp_theta_up1 = exp_theta_u * exp_theta[n];
        const T_partials_return exp_theta_vp1 = exp_theta_v * exp_theta[n];
        const T_partials_return exp_theta_uvsum = exp_theta_uv - exp_theta_up1 - exp_theta_vp1 + exp_theta[n];

        // Calculate the likelihood of Frank copula
        // if (include_summand<propto>::value)
        //   logp = 0;
        if (include_summand<propto, T_theta>::value)
            logp += log(theta_value[n] * (exp_theta[n]-1));


        if (include_summand<propto, T_u, T_v, T_theta>::value)
          logp += theta_value[n] * (u_dbl + v_dbl+1) -
                    2 * log(abs( exp_theta_uvsum )) ;
         // Calculate the derivative when the type is var (not double)
         if (!is_constant_struct<T_u>::value)
            ops_partials.edge1_.partials_[n] += theta_value[n] - 2 * theta_value[n] * ( exp_theta_uv - exp_theta_up1) / exp_theta_uvsum ;
         if (!is_constant_struct<T_v>::value)
            ops_partials.edge2_.partials_[n] += theta_value[n] - 2 * theta_value[n] * (exp_theta_uv - exp_theta_vp1) / exp_theta_uvsum ;
         if (!is_constant_struct<T_theta>::value)
            ops_partials.edge3_.partials_[n] += 1/ theta_value[n] + exp_theta[n] / (exp_theta[n] - 1) +
                                            (u_dbl + v_dbl+1) -
                                            2 * (exp_theta_uv * (u_dbl+v_dbl)- exp_theta_up1*(u_dbl+1) - exp_theta_vp1*(v_dbl+1) + exp_theta[n] )/
                                                exp_theta_uvsum;

      }
      return ops_partials.build(logp);
    }

    template <typename T_u, typename T_v, typename T_theta>
    inline
    typename stan::return_type<T_u, T_v, T_theta>::type
    bicop_r90_frank_log(const T_u& u, const T_v& v, const T_theta& theta) {
      return bicop_r90_frank_log<false>(u, v, theta);
    }

}
#endif // VIFCOPULA_DISTRIBUTION_BICOP_ROTATE_FRANK_LOG_HPP


