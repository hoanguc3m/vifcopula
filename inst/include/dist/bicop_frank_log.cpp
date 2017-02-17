#ifndef VIFCOPULA_DISTRIBUTION_BICOP_FRANK_LOG_CPP
#define VIFCOPULA_DISTRIBUTION_BICOP_FRANK_LOG_CPP

#include <Rcpp.h>
#include <stan/math.hpp>

// [[Rcpp::depends(StanHeaders)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(openmp)]]

namespace vifcopula {

using namespace Rcpp;
using namespace Eigen;
using namespace stan::math;
using namespace stan;

    /**
     * The log of the bivariate Frank copula density for the specified vector(s) u
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
    bicop_frank_log(const T_u& u, const T_v& v, const T_theta& theta) {
      static const char* function("vifcopula::bicop_frank_log");
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

      OperandsAndPartials<T_u, T_v, T_theta> operands_and_partials(u, v, theta);

      stan::VectorView<const T_u> u_vec(u);
      stan::VectorView<const T_v> v_vec(v);
      stan::VectorView<const T_theta> theta_vec(theta);

      size_t N = stan::max_size(u, v, theta);

      stan::VectorBuilder<true, T_partials_return, T_theta> theta_value(length(theta));
      stan::VectorBuilder<true, T_partials_return, T_theta> exp_neg_theta(length(theta));
      stan::VectorBuilder<true, T_partials_return, T_theta> om_exp_neg_theta(length(theta));
      stan::VectorBuilder<true, T_partials_return, T_theta> log_theta(length(theta));

      for (size_t i = 0; i < length(theta); i++) {
        theta_value[i] = value_of(theta_vec[i]);
        exp_neg_theta[i] = exp(- theta_value[i]);
        om_exp_neg_theta[i] = 1 - exp_neg_theta[i];
        log_theta[i] = log(theta_value[i]);
      }


      for (size_t n = 0; n < N; n++) {
        const T_partials_return u_dbl = value_of(u_vec[n]);
        const T_partials_return v_dbl = value_of(v_vec[n]);
        const T_partials_return exp_neg_theta_u = pow( exp_neg_theta[n], u_dbl);
        const T_partials_return exp_neg_theta_v = pow( exp_neg_theta[n], v_dbl);


        // Calculate the likelihood of Frank copula
        // if (include_summand<propto>::value)
        //   logp = 0;
        if (include_summand<propto, T_theta>::value)
            logp += log_theta[n] + log(om_exp_neg_theta[n]);


        if (include_summand<propto, T_u, T_v, T_theta>::value)
          logp -= theta_value[n] * (u_dbl + v_dbl) +
                    2 * log(om_exp_neg_theta[n] - (1 - exp_neg_theta_u) * (1 - exp_neg_theta_v) ) ;


//         // Calculate the derivative when the type is var (not double)
//         if (!is_constant_struct<T_u>::value)
//             operands_and_partials.d_x1[n] += - ( sq_theta[n] * inv_u_dbl[n] - theta_value[n] * inv_v_dbl[n] )/
//                                               (1 - sq_theta[n]) / pdf(s,inv_u_dbl[n])  ;
//         if (!is_constant_struct<T_v>::value)
//           operands_and_partials.d_x2[n] += - ( sq_theta[n] * inv_v_dbl[n] - theta_value[n] * inv_u_dbl[n] ) /
//                                               (1 - sq_theta[n]) / pdf(s,inv_v_dbl[n])  ;
//         if (!is_constant_struct<T_theta>::value)
//           operands_and_partials.d_x3[n] += theta_value[n] / (1 - sq_theta[n]) +
//                 ((1 + sq_theta[n]) * inv_u_dbl[n] * inv_v_dbl[n] - theta_value[n] * square(inv_u_dbl[n])
//                                                                - theta_value[n] * square(inv_v_dbl[n]) ) /
//                 square(1 - sq_theta[n]);
      }
      return operands_and_partials.value(logp);
    }

    template <typename T_u, typename T_v, typename T_theta>
    inline
    typename stan::return_type<T_u, T_v, T_theta>::type
    bicop_frank_log(const T_u& u, const T_v& v, const T_theta& theta) {
      return bicop_frank_log<false>(u, v, theta);
    }

}
#endif // VIFCOPULA_DISTRIBUTION_BICOP_FRANK_LOG_CPP

