#ifndef VIFCOPULA_DISTRIBUTION_BICOP_NORMAL_LOG_CPP
#define VIFCOPULA_DISTRIBUTION_BICOP_NORMAL_LOG_CPP

#include <stan/math.hpp>


namespace vifcopula {

using namespace stan::math;
using namespace stan;

    /**
     * The log of the bivariate normal copula density for the specified vector(s) u
     * and the specified vector(s) v given correlation(s) rho. u, v, or rho can
     * each be either a scalar or a vector. Any vector inputs
     * must be the same length.
     *
     * <p>The result log probability is defined to be the sum of the
     * log probabilities for each observation/observation/correlation triple.
     * @param u (Sequence of) scalar(s) in [0,1].
     * @param v (Sequence of) scalar(s) in [0,1].
     * @param rho (Sequence of) correlation parameters for the normal copula
     * @return The log of the product of the densities.
     * @tparam T_u Underlying type of scalar in sequence.
     * @tparam T_v Underlying type of scalar in sequence.
     * @tparam T_rho Type of correlation parameter.
     */

    template <bool propto, typename T_u, typename T_v, typename T_rho>
    typename stan::return_type<T_u, T_v, T_rho>::type
    bicop_normal_log(const T_u& u, const T_v& v, const T_rho& rho) {
      static const char* function("vifcopula::bicop_normal_log");
      typedef typename stan::partials_return_type<T_u, T_v, T_rho>::type
        T_partials_return;

      using std::log;
      using stan::is_constant_struct;
      using stan::math::square;
      using stan::math::log1m;
      using boost::math::normal;

      if (!(stan::length(u)
            && stan::length(v)
            && stan::length(rho)))
        return 0.0;

      T_partials_return logp(0.0);

//      check_bounded(function, "Random variable u", u,0,1);
//      check_bounded(function, "Random variable v", v,0,1);
//      check_bounded(function, "Correlation parameter", rho, -1,1);
      check_consistent_sizes(function,
                             "Random variable u", u,
                             "Random variable v", v,
                             "Scale parameter", rho);
      if (!include_summand<propto, T_u, T_v, T_rho>::value)
        return 0.0;

      OperandsAndPartials<T_u, T_v, T_rho> operands_and_partials(u, v, rho);

      stan::VectorView<const T_u> u_vec(u);
      stan::VectorView<const T_v> v_vec(v);
      stan::VectorView<const T_rho> rho_vec(rho);
      size_t N = stan::max_size(u, v, rho);

      stan::VectorBuilder<true, T_partials_return, T_rho> rho_value(length(rho));
      stan::VectorBuilder<true, T_partials_return, T_rho> sq_rho(length(rho));
      stan::VectorBuilder<include_summand<propto, T_rho>::value,
                    T_partials_return, T_rho> log_1mrhosq(length(rho));
      stan::VectorBuilder<true, T_partials_return, T_u> inv_u_dbl(length(u));
      stan::VectorBuilder<true, T_partials_return, T_v> inv_v_dbl(length(v));

      for (size_t i = 0; i < length(rho); i++) {
        rho_value[i] = value_of(rho_vec[i]);
        sq_rho[i] = square(rho_value[i]);
        if (include_summand<propto, T_rho>::value)
          log_1mrhosq[i] = log1m(square(rho_value[i]));
      }

      normal s;

      for (size_t n = 0; n < N; n++) {
        inv_u_dbl[n] = quantile(s,value_of(u_vec[n]));
        inv_v_dbl[n] = quantile(s,value_of(v_vec[n]));

        const T_partials_return iu_minus_iv_over_rho
          = (sq_rho[n] * ( square(inv_u_dbl[n]) + square(inv_v_dbl[n]) )
             - 2 * rho_value[n] * inv_u_dbl[n] * inv_v_dbl[n] ) / (1 - sq_rho[n]);

        // Calculate the likelihood of Gaussian copula
        static double NEGATIVE_HALF = - 0.5;
        // if (include_summand<propto>::value)
        //   logp += 0;
        if (include_summand<propto, T_rho>::value)
          logp += NEGATIVE_HALF * log_1mrhosq[n];
        if (include_summand<propto, T_u, T_v, T_rho>::value)
          logp += NEGATIVE_HALF * iu_minus_iv_over_rho;

        // Calculate the derivative when the type is var (not double)
        if (!is_constant_struct<T_u>::value)
            operands_and_partials.d_x1[n] += - ( sq_rho[n] * inv_u_dbl[n] - rho_value[n] * inv_v_dbl[n] )/
                                              (1 - sq_rho[n]) / pdf(s,inv_u_dbl[n])  ;
        if (!is_constant_struct<T_v>::value)
          operands_and_partials.d_x2[n] += - ( sq_rho[n] * inv_v_dbl[n] - rho_value[n] * inv_u_dbl[n] ) /
                                              (1 - sq_rho[n]) / pdf(s,inv_v_dbl[n])  ;
        if (!is_constant_struct<T_rho>::value)
          operands_and_partials.d_x3[n] += rho_value[n] / (1 - sq_rho[n]) +
                ((1 + sq_rho[n]) * inv_u_dbl[n] * inv_v_dbl[n] - rho_value[n] * square(inv_u_dbl[n])
                                                               - rho_value[n] * square(inv_v_dbl[n]) ) /
                square(1 - sq_rho[n]);
      }
      return operands_and_partials.value(logp);
    }

    template <typename T_u, typename T_v, typename T_rho>
    inline
    typename stan::return_type<T_u, T_v, T_rho>::type
    bicop_normal_log(const T_u& u, const T_v& v, const T_rho& rho) {
      return bicop_normal_log<false>(u, v, rho);
    }

}
#endif // VIFCOPULA_DISTRIBUTION_BICOP_NORMAL_LOG_CPP


