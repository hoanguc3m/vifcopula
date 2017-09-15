#ifndef VIFCOPULA_DISTRIBUTION_BICOP_CLAYTON_LOG_HPP
#define VIFCOPULA_DISTRIBUTION_BICOP_CLAYTON_LOG_HPP

#include <stan/math.hpp>

namespace vifcopula {

using namespace stan::math;
using namespace stan;

    /**
     * The log of the bivariate Clayton copula density for the specified vector(s) u
     * and the specified vector(s) v given degree(s) of dependence theta.
     *  u, v, or theta can each be either a scalar or a vector. Any vector inputs
     * must be the same length.
     *
     * <p>The result log probability is defined to be the sum of the
     * log probabilities for each observation/observation/dependence.
     * @param u (Sequence of) scalar(s) in [0,1].
     * @param v (Sequence of) scalar(s) in [0,1].
     * @param theta (Sequence of) dependence parameters for the Clayton copula
     * @return The log of the product of the densities.
     * @tparam T_u Underlying type of scalar in sequence.
     * @tparam T_v Underlying type of scalar in sequence.
     * @tparam T_theta Type of dependence parameter.
     */

    template <bool propto, typename T_u, typename T_v, typename T_theta>
    typename stan::return_type<T_u, T_v, T_theta>::type
    bicop_clayton_log(const T_u& u, const T_v& v, const T_theta& theta) {
      static const char* function("vifcopula::bicop_clayton_log");
      typedef typename stan::partials_return_type<T_u, T_v, T_theta>::type
        T_partials_return;

      using std::log;

      if (!(stan::length(u)
            && stan::length(v)
            && stan::length(theta)))
        return 0.0;

      T_partials_return logp(0.0);

//      check_bounded(function, "Random variable u", u,0,1);
//      check_bounded(function, "Random variable v", v,0,1);
//      check_positive_finite(function, "dependence parameter", theta);
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
      stan::VectorBuilder<true, T_partials_return, T_theta> theta_p1(length(theta));
      stan::VectorBuilder<true, T_partials_return, T_theta> inv_theta_p2(length(theta));

      for (size_t i = 0; i < length(theta); i++) {
        theta_value[i] = value_of(theta_vec[i]);
        inv_theta[i] = 1 / theta_value[i];
        theta_p1[i] = 1 + theta_value[i];
        inv_theta_p2[i] = inv_theta[i] + 2;
      }


      for (size_t n = 0; n < N; n++) {
        T_partials_return u_dbl = value_of(u_vec[n]);        // No const
        T_partials_return v_dbl = value_of(v_vec[n]);        // No const
        
        if (u_dbl < 1e-10) { u_dbl = 1e-10;}
        if (v_dbl < 1e-10) { v_dbl = 1e-10;}
        if (u_dbl > 1-1e-10) { u_dbl = 1-1e-10;}
        if (v_dbl > 1-1e-10) { v_dbl = 1-1e-10;}
        
        const T_partials_return log_uv = log(u_dbl) + log(v_dbl);
        T_partials_return A_u_v_theta
                    = pow(u_dbl,-theta_value[n]) + pow(v_dbl,-theta_value[n]) -1 ; // No const
        // if A_u_v_theta = Inf
        if (A_u_v_theta >  std::numeric_limits<double>::max() ){
            A_u_v_theta = std::numeric_limits<double>::max();
        }


        const T_partials_return log_A = log(A_u_v_theta);

        // Calculate the likelihood of clayton copula
        // if (include_summand<propto>::value)
        //   logp = 0;

        if (include_summand<propto, T_theta>::value)
          logp += log(theta_p1[n]);

        if (include_summand<propto, T_u, T_v, T_theta>::value)
          logp -= theta_p1[n] * log_uv + inv_theta_p2[n] * log_A ;


         // Calculate the derivative when the type is var (not double)
         if (!is_constant_struct<T_u>::value)
             ops_partials.edge1_.partials_[n] +=  - theta_p1[n] / u_dbl +  theta_value[n] * inv_theta_p2[n] /u_dbl/( 1 + pow(v_dbl/u_dbl,-theta_value[n]) - pow(u_dbl,theta_value[n])  );
         if (!is_constant_struct<T_v>::value)
            ops_partials.edge2_.partials_[n] +=  - theta_p1[n] / v_dbl +  theta_value[n] * inv_theta_p2[n] /v_dbl/( 1 + pow(u_dbl/v_dbl,-theta_value[n]) - pow(v_dbl,theta_value[n])  );
         

         if (!is_constant_struct<T_theta>::value){
            T_partials_return A_u_v_theta_log_div_A
                    = (log(u_dbl) + pow(v_dbl/u_dbl,-theta_value[n]) * log(v_dbl) ) / ( 1 + pow(v_dbl/u_dbl,-theta_value[n]) - pow(u_dbl,theta_value[n])  );
            if (A_u_v_theta_log_div_A >  std::numeric_limits<double>::max() ){
                A_u_v_theta_log_div_A = std::max(log(v_dbl), log(v_dbl));
            } // check maximum double 
            ops_partials.edge3_.partials_[n] += 1/theta_p1[n] - log_uv + log_A/ square(theta_value[n]) +
                                                inv_theta_p2[n] * A_u_v_theta_log_div_A ;

         }


      }
      return ops_partials.build(logp);
    }

    template <typename T_u, typename T_v, typename T_theta>
    inline
    typename stan::return_type<T_u, T_v, T_theta>::type
    bicop_clayton_log(const T_u& u, const T_v& v, const T_theta& theta) {
      return bicop_clayton_log<false>(u, v, theta);
    }

}
#endif // VIFCOPULA_DISTRIBUTION_BICOP_CLAYTON_LOG_HPP


