#ifndef VIFCOPULA_DISTRIBUTION_BICOP_JOE_LOG_HPP
#define VIFCOPULA_DISTRIBUTION_BICOP_JOE_LOG_HPP

#include <stan/math.hpp>

namespace vifcopula {

using namespace stan::math;
using namespace stan;

    /**
     * The log of the bivariate Joe copula density for the specified vector(s) u
     * and the specified vector(s) v given degree(s) of dependence theta.
     *  u, v, or theta can each be either a scalar or a vector. Any vector inputs
     * must be the same length.
     *
     * <p>The result log probability is defined to be the sum of the
     * log probabilities for each observation/observation/dependence.
     * @param u (Sequence of) scalar(s) in [0,1].
     * @param v (Sequence of) scalar(s) in [0,1].
     * @param theta (Sequence of) dependence parameters for the Joe copula
     * @return The log of the product of the densities.
     * @tparam T_u Underlying type of scalar in sequence.
     * @tparam T_v Underlying type of scalar in sequence.
     * @tparam T_theta Type of dependence parameter.
     */

    template <bool propto, typename T_u, typename T_v, typename T_theta>
    typename stan::return_type<T_u, T_v, T_theta>::type
    bicop_joe_log(const T_u& u, const T_v& v, const T_theta& theta) {
      static const char* function("vifcopula::bicop_joe_log");
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
//      check_positive(function, "dependence parameter", theta - 1); // theta > 1
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
      stan::VectorBuilder<true, T_partials_return, T_theta> inv_theta(length(theta));
      stan::VectorBuilder<true, T_partials_return, T_theta> thetam1(length(theta));

      for (size_t i = 0; i < length(theta); i++) {
        theta_value[i] = value_of(theta_vec[i]);
        inv_theta[i] = 1 / theta_value[i];
        thetam1[i] = theta_value[i] -1;
      }


      for (size_t n = 0; n < N; n++) {
        const T_partials_return u_dbl = value_of(u_vec[n]);
        const T_partials_return v_dbl = value_of(v_vec[n]);

        const T_partials_return om_u = 1 - u_dbl;
        const T_partials_return om_v = 1 - v_dbl;
        const T_partials_return t_u = pow(om_u, theta_value[n]);
        const T_partials_return t_v = pow(om_v, theta_value[n]);
        const T_partials_return t_uv = t_u + t_v - t_u * t_v;

        const T_partials_return dtu_theta = t_u * log(om_u);
        const T_partials_return dtv_theta = t_v * log(om_v);
        const T_partials_return dtuv_sum = dtu_theta + dtv_theta - dtu_theta * t_v - dtv_theta*t_u ;

        // Calculate the likelihood of Joe copula
        // if (include_summand<propto>::value)
        //   logp = 0;

        if (include_summand<propto, T_u, T_v, T_theta>::value)
          logp += log(t_uv) * (inv_theta[n]-2) +
                    log(thetam1[n] + t_uv) + thetam1[n] * log(om_u) + thetam1[n] * log(om_v);


         // Calculate the derivative when the type is var (not double)
         if (!is_constant_struct<T_u>::value)
             ops_partials.edge1_.partials_[n] += (inv_theta[n] - 2) * (t_v-1) * theta_value[n] * t_u / (1-u_dbl)/t_uv +
                                                (t_v-1) * theta_value[n] * t_u / (1-u_dbl)/(t_uv + thetam1[n]) - thetam1[n]/om_u;
         if (!is_constant_struct<T_v>::value)
           ops_partials.edge2_.partials_[n] +=  (inv_theta[n] - 2) * (t_u-1) * theta_value[n] * t_v / (1-v_dbl)/t_uv +
                                                (t_u-1) * theta_value[n] * t_v / (1-v_dbl)/(t_uv + thetam1[n]) - thetam1[n]/om_v;
         if (!is_constant_struct<T_theta>::value)
           ops_partials.edge3_.partials_[n] +=  - log(t_uv)/square(theta_value[n]) +
                                        (inv_theta[n] - 2)*dtuv_sum/t_uv +
                                        log(om_u) + log(om_v) + (1 + dtuv_sum)/(thetam1[n] + t_uv) ;
      }
      return ops_partials.build(logp);
    }

    template <typename T_u, typename T_v, typename T_theta>
    inline
    typename stan::return_type<T_u, T_v, T_theta>::type
    bicop_joe_log(const T_u& u, const T_v& v, const T_theta& theta) {
      return bicop_joe_log<false>(u, v, theta);
    }

}
#endif // VIFCOPULA_DISTRIBUTION_BICOP_JOE_LOG_HPP


