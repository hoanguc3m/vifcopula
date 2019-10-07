#ifndef VIFCOPULA_DISTRIBUTION_BICOP_SURVIVAL_BB1_LOG_HPP
#define VIFCOPULA_DISTRIBUTION_BICOP_SURVIVAL_BB1_LOG_HPP

#include <stan/math.hpp>
#include <iostream>

namespace vifcopula {

using namespace stan::math;
using namespace stan;

    /**
     * The log of the bivariate SURVIVAL BB1 copula density for the specified vector(s) u
     * and the specified vector(s) v given degree(s) of dependence theta.
     *  u, v, or theta can each be either a scalar or a vector. Any vector inputs
     * must be the same length.
     *
     * <p>The result log probability is defined to be the sum of the
     * log probabilities for each observation/observation/dependence.
     * @param u (Sequence of) scalar(s) in [0,1].
     * @param v (Sequence of) scalar(s) in [0,1].
     * @param theta (Sequence of) dependence parameters for the BB1 copula
     * @return The log of the product of the densities.
     * @tparam T_u Underlying type of scalar in sequence.
     * @tparam T_v Underlying type of scalar in sequence.
     * @tparam T_theta Type of dependence parameter.
     */

    template <bool propto, typename T_u, typename T_v, typename T_theta, typename T_delta>
    typename stan::return_type<T_u, T_v, T_theta, T_delta>::type
    bicop_survival_BB1_log(const T_u& u, const T_v& v, const T_theta& theta, const T_delta& delta) {
      static const char* function("vifcopula::bicop_survival_BB1_log");
      typedef typename stan::partials_return_type<T_u, T_v, T_theta, T_delta>::type
        T_partials_return;

      using std::log;

      if (!(stan::length(u)
            && stan::length(v)
            && stan::length(theta)
            && stan::length(delta)))
        return 0.0;

      T_partials_return logp(0.0);

//      check_bounded(function, "Random variable u", u,0,1);
//      check_bounded(function, "Random variable v", v,0,1);
//      check_positive_finite(function, "dependence parameter", theta);
//      check_bounded(function, "dependence parameter", delta,0,1);
      check_consistent_sizes(function,
                             "Random variable u", u,
                             "Random variable v", v,
                             "theta parameter", theta,
                             "theta parameter", delta);
      if (!include_summand<propto, T_u, T_v, T_theta, T_delta>::value)
        return 0.0;

      operands_and_partials<T_u, T_v, T_theta, T_delta> ops_partials(u, v, theta, delta);

      scalar_seq_view<T_u> u_vec(u);
      scalar_seq_view<T_v> v_vec(v);
      scalar_seq_view<T_theta> theta_vec(theta);
      scalar_seq_view<T_delta> delta_vec(delta);

      size_t N = stan::max_size(u, v, theta, delta);

      stan::VectorBuilder<true, T_partials_return, T_theta> theta_value(length(theta));
      stan::VectorBuilder<true, T_partials_return, T_delta> delta_value(length(delta));

      stan::VectorBuilder<true, T_partials_return, T_theta> inv_theta(length(theta));
      stan::VectorBuilder<true, T_partials_return, T_delta> inv_delta(length(delta));


      for (size_t i = 0; i < length(theta); i++) {
        theta_value[i] = value_of(theta_vec[i]);
        delta_value[i] = value_of(delta_vec[i]);

        inv_theta[i] = 1 / theta_value[i];
        inv_delta[i] = 1 / delta_value[i];
      }


      for (size_t n = 0; n < N; n++) {
        T_partials_return u_dbl = 1. - value_of(u_vec[n]);        // survial
        T_partials_return v_dbl = 1. - value_of(v_vec[n]);        // survial

        if (u_dbl < 1e-10) { u_dbl = 1e-10;}
        if (v_dbl < 1e-10) { v_dbl = 1e-10;}
        if (u_dbl > 1-1e-10) { u_dbl = 1-1e-10;}
        if (v_dbl > 1-1e-10) { v_dbl = 1-1e-10;}

        T_partials_return x_dbl = pow(pow(u_dbl,-theta_value[n]) - 1,delta_value[n]);
        T_partials_return y_dbl = pow(pow(v_dbl,-theta_value[n]) - 1,delta_value[n]);

        T_partials_return logxy = log(x_dbl) + log(y_dbl);
        T_partials_return xpy = x_dbl + y_dbl;

        T_partials_return xpyd = pow(xpy, inv_delta[n]);
        T_partials_return thde = theta_value[n] * delta_value[n];

        T_partials_return A_xpydeth = theta_value[n] * (delta_value[n]-1) + (thde + 1)* pow(xpy, inv_delta[n]) ;
        T_partials_return log_c = -(inv_theta[n]+2) * log(1 + xpyd) +\
                                   (inv_delta[n]-2) * log(xpy) + \
                                    log(A_xpydeth) +\
                                    (1 - inv_delta[n]) * logxy -\
                                    (theta_value[n] + 1) * (log(u_dbl) + log(v_dbl)) ;

        if (include_summand<propto, T_u, T_v, T_theta, T_delta>::value)
          logp += log_c;



         // Calculate the derivative when the type is var (not double)
         // survial  - negative the deriv
         if (!is_constant_struct<T_u>::value){
           T_partials_return dx_du = - thde * pow(pow(u_dbl,-theta_value[n]) - 1.,delta_value[n]-1.) * pow(u_dbl,-theta_value[n]-1.);
           T_partials_return xpydm1 =  xpyd / xpy;
           ops_partials.edge1_.partials_[n] -=  -(theta_value[n] + 1.) / u_dbl + \
             dx_du * ( - (1./thde + 2.*inv_delta[n]) * xpydm1 / (1. + xpyd) + \
                       (inv_delta[n] - 2) /  xpy  + \
                       (theta_value[n] + inv_delta[n]) * xpydm1 / A_xpydeth +\
                       (1 - inv_delta[n]) / x_dbl ) ;

         }


         if (!is_constant_struct<T_v>::value){
           T_partials_return dy_dv = - thde * pow(pow(v_dbl,-theta_value[n]) - 1.,delta_value[n]-1.) * pow(v_dbl,-theta_value[n]-1.);
           T_partials_return xpydm1 =  xpyd / xpy;
           ops_partials.edge2_.partials_[n] -=  -(theta_value[n] + 1.) / v_dbl + \
             dy_dv * ( - (1./thde + 2.*inv_delta[n]) * xpydm1 / (1. + xpyd) + \
                        (inv_delta[n] - 2.) /  xpy  + \
                        (theta_value[n] + inv_delta[n]) * xpydm1 / A_xpydeth + \
                        (1 - inv_delta[n]) / y_dbl ) ;

        }

         if (!is_constant_struct<T_theta>::value){
           double eps = 0.001;
           T_partials_return theta_new = theta_value[n] + eps;
           T_partials_return delta_new = delta_value[n];

           T_partials_return inv_theta_new = 1 / theta_new;
           T_partials_return inv_delta_new = 1 / delta_new;

           T_partials_return x_new = pow(pow(u_dbl,-theta_new) - 1,delta_new);
           T_partials_return y_new = pow(pow(v_dbl,-theta_new) - 1,delta_new);

           T_partials_return logxy_new = log(x_new) + log(y_new);
           T_partials_return xpy_new = x_new + y_new;

           T_partials_return xpyd_new = pow(xpy_new, inv_delta_new);
           T_partials_return thde_new = theta_new * delta_new;

           T_partials_return A_xpydeth_new = theta_new * (delta_new-1) + (thde_new + 1)* pow(xpy_new, inv_delta_new) ;
           T_partials_return log_c_new = -(inv_theta_new + 2) * log(1 + xpyd_new) +\
             (inv_delta_new-2) * log(xpy_new) +                              \
             log(A_xpydeth_new) +                                           \
             (1 - inv_delta_new) * logxy_new -                               \
             (theta_new + 1) * (log(u_dbl) + log(v_dbl)) ;
            ops_partials.edge3_.partials_[n] += (log_c_new - log_c)/eps;
         }

         if (!is_constant_struct<T_delta>::value){
           double eps = 0.001;
           T_partials_return theta_new = theta_value[n];
           T_partials_return delta_new = delta_value[n] + eps;

           T_partials_return inv_theta_new = 1 / theta_new;
           T_partials_return inv_delta_new = 1 / delta_new;

           T_partials_return x_new = pow(pow(u_dbl,-theta_new) - 1,delta_new);
           T_partials_return y_new = pow(pow(v_dbl,-theta_new) - 1,delta_new);

           T_partials_return logxy_new = log(x_new) + log(y_new);
           T_partials_return xpy_new = x_new + y_new;

           T_partials_return xpyd_new = pow(xpy_new, inv_delta_new);
           T_partials_return thde_new = theta_new * delta_new;

           T_partials_return A_xpydeth_new = theta_new * (delta_new-1) + (thde_new + 1)* pow(xpy_new, inv_delta_new) ;
           T_partials_return log_c_new = -(inv_theta_new + 2) * log(1 + xpyd_new) +\
             (inv_delta_new-2) * log(xpy_new) +                              \
             log(A_xpydeth_new) +                                           \
             (1 - inv_delta_new) * logxy_new -                               \
             (theta_new + 1) * (log(u_dbl) + log(v_dbl)) ;
           ops_partials.edge4_.partials_[n] += (log_c_new - log_c)/eps;

         }


      }
      return ops_partials.build(logp);
    }

    template <typename T_u, typename T_v, typename T_theta, typename T_delta>
    inline
    typename stan::return_type<T_u, T_v, T_theta, T_delta>::type
    bicop_survival_BB1_log(const T_u& u, const T_v& v, const T_theta& theta, const T_delta& delta) {
      return bicop_survival_BB1_log<false>(u, v, theta, delta);
    }

}
#endif // VIFCOPULA_DISTRIBUTION_BICOP_SURVIVAL_BB1_LOG_HPP


