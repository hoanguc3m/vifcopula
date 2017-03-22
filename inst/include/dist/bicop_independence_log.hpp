#ifndef VIFCOPULA_DISTRIBUTION_BICOP_INDEPENDENCE_LOG_HPP
#define VIFCOPULA_DISTRIBUTION_BICOP_INDEPENDENCE_LOG_HPP

#include <stan/math.hpp>


namespace vifcopula {

using namespace stan::math;
using namespace stan;

    /**
     * The log of the bivariate independence copula density for the specified
     * vector(s) u and the specified vector(s) v. u or v can
     * each be either a scalar or a vector. Any vector inputs
     * must be the same length.
     *
     * <p>The result log probability is defined to be the sum of the
     * log probabilities for each observation/observation.
     * @param u (Sequence of) scalar(s) in [0,1].
     * @param v (Sequence of) scalar(s) in [0,1].
     * @return The log of the product of the densities.
     * @tparam T_u Underlying type of scalar in sequence.
     * @tparam T_v Underlying type of scalar in sequence.
     */

    template <bool propto, typename T_u, typename T_v>
    typename stan::return_type<T_u, T_v>::type
    bicop_independence_log(const T_u& u, const T_v& v) {
      static const char* function("vifcopula::bicop_independence_log");
      typedef typename stan::partials_return_type<T_u, T_v>::type
        T_partials_return;

      using std::log;
      using stan::is_constant_struct;

      if (!(stan::length(u)
            && stan::length(v)))
        return 0.0;

      T_partials_return logp(0.0);

//      check_bounded(function, "Random variable u", u,0,1);
//      check_bounded(function, "Random variable v", v,0,1);
      check_consistent_sizes(function,
                             "Random variable u", u,
                             "Random variable v", v);
      if (!include_summand<propto, T_u, T_v>::value)
        return 0.0;

      OperandsAndPartials<T_u, T_v> operands_and_partials(u, v);

      stan::VectorView<const T_u> u_vec(u);
      stan::VectorView<const T_v> v_vec(v);
      size_t N = stan::max_size(u, v);


      for (size_t n = 0; n < N; n++) {

        // Calculate the likelihood of independence copula
        // log c(u,v) = 0


        // Calculate the derivative
        if (!is_constant_struct<T_u>::value)
            operands_and_partials.d_x1[n] += 1;
        if (!is_constant_struct<T_v>::value)
          operands_and_partials.d_x2[n] += 1;
      }
      return operands_and_partials.value(logp);
    }

    template <typename T_u, typename T_v>
    inline
    typename stan::return_type<T_u, T_v>::type
        bicop_independence_log(const T_u& u, const T_v& v) {
      return bicop_independence_log<false>(u, v);
    }

}
#endif // VIFCOPULA_DISTRIBUTION_BICOP_INDEPENDENCE_LOG_HPP


