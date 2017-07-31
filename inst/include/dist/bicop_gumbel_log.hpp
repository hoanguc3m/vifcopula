#ifndef VIFCOPULA_DISTRIBUTION_BICOP_GUMBEL_LOG_HPP
#define VIFCOPULA_DISTRIBUTION_BICOP_GUMBEL_LOG_HPP

#include <stan/math.hpp>

namespace vifcopula {

using namespace stan::math;
using namespace stan;

    /**
     * The log of the bivariate Gumbel copula density for the specified vector(s) u
     * and the specified vector(s) v given degree(s) of dependence theta.
     *  u, v, or theta can each be either a scalar or a vector. Any vector inputs
     * must be the same length.
     *
     * <p>The result log probability is defined to be the sum of the
     * log probabilities for each observation/observation/dependence.
     * @param u (Sequence of) scalar(s) in [0,1].
     * @param v (Sequence of) scalar(s) in [0,1].
     * @param theta (Sequence of) dependence parameters for the Gumbel copula
     * @return The log of the product of the densities.
     * @tparam T_u Underlying type of scalar in sequence.
     * @tparam T_v Underlying type of scalar in sequence.
     * @tparam T_theta Type of dependence parameter.
     */

    template <bool propto, typename T_u, typename T_v, typename T_theta>
    typename stan::return_type<T_u, T_v, T_theta>::type
    bicop_gumbel_log(const T_u& u, const T_v& v, const T_theta& theta) {
      static const char* function("vifcopula::bicop_gumbel_log");
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
//      check_positive_finite(function, "dependence parameter", theta - 1); // theta >1
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
      stan::VectorBuilder<true, T_partials_return, T_theta> theta_m1(length(theta));
      stan::VectorBuilder<true, T_partials_return, T_theta> inv_theta_t2m2(length(theta));

      for (size_t i = 0; i < length(theta); i++) {
        theta_value[i] = value_of(theta_vec[i]);
        inv_theta[i] = 1 / theta_value[i];
        theta_m1[i] = theta_value[i] -1;
        inv_theta_t2m2[i] = inv_theta[i] * 2 - 2;
      }


      for (size_t n = 0; n < N; n++) {
        const T_partials_return u_dbl = value_of(u_vec[n]);
        const T_partials_return v_dbl = value_of(v_vec[n]);
        const T_partials_return t_u = pow(-log(u_dbl), theta_value[n]);
        const T_partials_return t_v = pow(-log(v_dbl), theta_value[n]);
        const T_partials_return t_uv = t_u + t_v;
        const T_partials_return t_u_div_t_uv = t_u / t_uv;
        const T_partials_return t_v_div_t_uv = t_v / t_uv;
        const T_partials_return inv_u_logu = 1/(u_dbl * log(u_dbl)) ;
        const T_partials_return inv_v_logv = 1/(v_dbl * log(v_dbl)) ;

        const T_partials_return log_C_u_v_theta
                    = - pow(t_uv,inv_theta[n]);
        const T_partials_return op_theta_tuv = 1 + theta_m1[n] * pow(t_uv, - inv_theta[n] );

        // Calculate the likelihood of Gumbel copula
        // if (include_summand<propto>::value)
        //   logp = 0;

        if (include_summand<propto, T_u, T_v, T_theta>::value)
          logp += log_C_u_v_theta - log(u_dbl * v_dbl) +
                    log(t_uv) * inv_theta_t2m2[n] + (theta_m1[n]) * log( log(u_dbl) * log(v_dbl) ) +
                    log( op_theta_tuv );


         // Calculate the derivative when the type is var (not double)
         if (!is_constant_struct<T_u>::value)
             ops_partials.edge1_.partials_[n] +=  log_C_u_v_theta * t_u_div_t_uv * inv_u_logu - 1/u_dbl +
                                            inv_theta_t2m2[n] * theta_value[n] * t_u_div_t_uv * inv_u_logu +
                                            theta_m1[n] * inv_u_logu +
                                            theta_m1[n] * t_u_div_t_uv * inv_u_logu / log_C_u_v_theta / op_theta_tuv ;
         if (!is_constant_struct<T_v>::value)
           ops_partials.edge2_.partials_[n] +=  log_C_u_v_theta * t_v_div_t_uv * inv_v_logv - 1/v_dbl +
                                            inv_theta_t2m2[n] * theta_value[n] * t_v_div_t_uv * inv_v_logv +
                                            theta_m1[n] * inv_v_logv +
                                            theta_m1[n] * t_v_div_t_uv * inv_v_logv / log_C_u_v_theta / op_theta_tuv ;

         if (!is_constant_struct<T_theta>::value){
            const T_partials_return log_u           = log(u_dbl);
			const T_partials_return log_v           = log(v_dbl);
			const T_partials_return C_u_v_theta     = exp(log_C_u_v_theta);
			const T_partials_return logtuv_div_theta_sq = log(t_uv)/(theta_value[n]*theta_value[n]);
            const T_partials_return t_u_log_t_v_log = t_u*log(-log_u)+t_v*log(-log_v);

			const T_partials_return log_tuv_tuv     = -logtuv_div_theta_sq+inv_theta[n]*t_u_log_t_v_log/t_uv;

			const T_partials_return t_uv_pow_thetat2m2 = pow(t_uv,inv_theta_t2m2[n]);

			const T_partials_return C_u_v_t_uv      = C_u_v_theta*t_uv_pow_thetat2m2;
			const T_partials_return log_u_log_v     = log(u_dbl)*log(v_dbl);

			const T_partials_return log_umv_thetam1 = pow(log_u_log_v,theta_m1[n]);
			const T_partials_return t_uv_pow_minv_theta = pow(t_uv,-inv_theta[n]);
			const T_partials_return o_div_uv        = 1/(u_dbl*v_dbl);
			const T_partials_return log_umv_opo = log_umv_thetam1*op_theta_tuv*o_div_uv;

			const T_partials_return C_uvt_log_umv   = C_u_v_t_uv*log_umv_thetam1;

           ops_partials.edge3_.partials_[n] += (log_C_u_v_theta*log_tuv_tuv*C_u_v_t_uv*log_umv_opo
                                                + C_u_v_t_uv*(-2.0*logtuv_div_theta_sq +inv_theta_t2m2[n]*t_u_log_t_v_log/t_uv)*log_umv_opo
                                                +C_uvt_log_umv*log(log_u_log_v)*op_theta_tuv*o_div_uv
                                                +C_uvt_log_umv*(t_uv_pow_minv_theta-(op_theta_tuv-1)*log_tuv_tuv)*o_div_uv) /
                                                C_u_v_theta/t_uv_pow_thetat2m2/log_umv_thetam1/op_theta_tuv*u_dbl*v_dbl;
         }
      }
      return ops_partials.build(logp);
    }

    template <typename T_u, typename T_v, typename T_theta>
    inline
    typename stan::return_type<T_u, T_v, T_theta>::type
    bicop_gumbel_log(const T_u& u, const T_v& v, const T_theta& theta) {
      return bicop_gumbel_log<false>(u, v, theta);
    }

}
#endif // VIFCOPULA_DISTRIBUTION_BICOP_GUMBEL_LOG_HPP


