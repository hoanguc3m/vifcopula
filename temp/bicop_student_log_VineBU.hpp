#ifndef VIFCOPULA_DISTRIBUTION_BICOP_STUDENT_LOG_HPP
#define VIFCOPULA_DISTRIBUTION_BICOP_STUDENT_LOG_HPP

#include <stan/math.hpp>
#include <iostream>

void (*difflPDF_nu_tCopula_new) (double* u, double* v, int* n, double* param, int* copula, double* out);    // deriv par2

namespace vifcopula {

using namespace stan::math;
using namespace stan;

    /**
     * The log of the bivariate Student copula density for the specified vector(s) u
     * and the specified vector(s) v given correlation(s) rho and degree of freedom nu.
     *  u, v, or rho, nu can each be either a scalar or a vector. Any vector inputs
     * must be the same length.
     *
     * <p>The result log probability is defined to be the sum of the
     * log probabilities for each observation/observation/correlation/df.
     * @param u (Sequence of) scalar(s) in [0,1].
     * @param v (Sequence of) scalar(s) in [0,1].
     * @param rho (Sequence of) correlation parameters for the Student copula
     * @param nu (Sequence of) freedom parameters for the Student copula
     * @return The log of the product of the densities.
     * @tparam T_u Underlying type of scalar in sequence.
     * @tparam T_v Underlying type of scalar in sequence.
     * @tparam T_rho Type of correlation parameter.
     * @tparam T_nu Type of correlation parameter.
     */

    template <bool propto, typename T_u, typename T_v, typename T_rho, typename T_nu>
    typename stan::return_type<T_u, T_v, T_rho, T_nu>::type
    bicop_student_log(const T_u& u, const T_v& v, const T_rho& rho, const T_nu& nu ) {
      static const char* function("vifcopula::bicop_student_log");
      typedef typename stan::partials_return_type<T_u, T_v, T_rho, T_nu>::type
        T_partials_return;

      using std::log;
      using stan::math::lgamma;
      using stan::math::digamma;
      using stan::math::lbeta;
      using stan::math::grad_reg_inc_beta;
      using boost::math::students_t;


      if (!(stan::length(u)
            && stan::length(v)
            && stan::length(rho)
            && stan::length(nu)))
        return 0.0;

      T_partials_return logp(0.0);

//      check_bounded(function, "Random variable u", u,0,1);
//      check_bounded(function, "Random variable v", v,0,1);
//      check_bounded(function, "Correlation parameter", rho, -1,1);
//      check_positive_finite(function, "Degree of freedom parameter", nu);
      check_consistent_sizes(function,
                             "Random variable u", u,
                             "Random variable v", v,
                             "Scale parameter", rho,
                             "Degree of freedom parameter", nu);
      if (!include_summand<propto, T_u, T_v, T_rho, T_nu>::value)
        return 0.0;

      operands_and_partials<T_u, T_v, T_rho, T_nu> ops_partials(u, v, rho, nu);

      scalar_seq_view<T_u> u_vec(u);
      scalar_seq_view<T_v> v_vec(v);
      scalar_seq_view<T_rho> rho_vec(rho);
      scalar_seq_view<T_nu> nu_vec(nu);
      size_t N = stan::max_size(u, v, rho, nu);

      stan::VectorBuilder<true, T_partials_return, T_rho> rho_value(length(rho));
      stan::VectorBuilder<true, T_partials_return, T_rho> sq_rho(length(rho));
      stan::VectorBuilder<include_summand<propto, T_rho>::value,
                    T_partials_return, T_rho> log_1mrhosq(length(rho));

      stan::VectorBuilder<true, T_partials_return, T_nu> nu_value(length(nu));
      stan::VectorBuilder<true, T_partials_return, T_nu> nud2(length(nu));
      stan::VectorBuilder<true, T_partials_return, T_nu> nud2ph(length(nu));
      stan::VectorBuilder<true, T_partials_return, T_nu> nud2m1(length(nu));
      stan::VectorBuilder<true, T_partials_return, T_nu> nud2p1(length(nu));

      stan::VectorBuilder<true, T_partials_return, T_u> inv_u_dbl(length(u));
      stan::VectorBuilder<true, T_partials_return, T_v> inv_v_dbl(length(v));

      for (size_t i = 0; i < length(rho); i++) {
        rho_value[i] = value_of(rho_vec[i]);
        nu_value[i] = value_of(nu_vec[i]);
        sq_rho[i] = square(rho_value[i]);
        if (include_summand<propto, T_rho>::value)
          log_1mrhosq[i] = log1m(square(rho_value[i]));
        nud2[i] = nu_value[i] * 0.5;
        nud2ph[i] = nud2[i] + 0.5;
        nud2m1[i] = nud2[i] - 1.0;
        nud2p1[i] = nud2[i] + 1.0;
      }


      for (size_t n = 0; n < N; n++) {
        students_t s(nu_value[n]);
        inv_u_dbl[n] = quantile(s,value_of(u_vec[n]));
        inv_v_dbl[n] = quantile(s,value_of(v_vec[n]));

        const T_partials_return M_nu_rho
                    = nu_value[n] * (1 - sq_rho[n]) +
                        square(inv_u_dbl[n]) + square(inv_v_dbl[n]) -
                        2 * rho_value[n] * inv_u_dbl[n] * inv_v_dbl[n]  ;

        const T_partials_return logM = log ( M_nu_rho ) ;

        // Calculate the likelihood of Student copula
        static double NEGATIVE_HALF = - 0.5;
        static double POSITIVE_HALF = 0.5;

        if (include_summand<propto>::value)
          logp -= log(2.0);

        if (include_summand<propto, T_nu>::value)
          logp += 2.0 * lgamma(nud2[n]) - 2.0 * lgamma(nud2ph[n]) - nud2m1[n] * log(nu_value[n]);

        if (include_summand<propto, T_rho, T_nu>::value)
          logp += nud2ph[n] * log_1mrhosq[n];

        if (include_summand<propto, T_u, T_v, T_rho, T_nu>::value)
          logp += nud2ph[n] * ( log(nu_value[n] + square(inv_u_dbl[n])) +
                                    log(nu_value[n] + square(inv_v_dbl[n]))) -
                                    nud2p1[n] *   logM              ;

         // Calculate the derivative when the type is var (not double)
         if (!is_constant_all<T_u>::value)
            ops_partials.edge1_.partials_[n] += ((nu_value[n] + 1) *  inv_u_dbl[n] / (  nu_value[n] + square(inv_u_dbl[n]) )
                                                - (nu_value[n] + 2) * (inv_u_dbl[n] - rho_value[n] * inv_v_dbl[n] / M_nu_rho ) ) /
                                                pdf(s,inv_u_dbl[n]) ;
         if (!is_constant_all<T_v>::value)
            ops_partials.edge2_.partials_[n] += ((nu_value[n] + 1) *  inv_v_dbl[n] / (  nu_value[n] + square(inv_v_dbl[n]) )
                                            - (nu_value[n] + 2) * (inv_v_dbl[n] - rho_value[n] * inv_u_dbl[n] ) / M_nu_rho ) /
                                            pdf(s,inv_v_dbl[n]);

         if (!is_constant_all<T_rho>::value)
            ops_partials.edge3_.partials_[n] += - (nu_value[n] + 1) * rho_value[n] / (1 - sq_rho[n]) +
                                               (nu_value[n] + 2) * (nu_value[n] * rho_value[n] + inv_u_dbl[n] * inv_v_dbl[n] ) / M_nu_rho  ;

         if (!is_constant_all<T_nu>::value) {
             double u_double = value_of(u_vec[n]);
             double v_double = value_of(v_vec[n]);
             double rho_double = value_of(rho_vec[n]);
             double nu_double = value_of(nu_vec[n]);
             int nn = 1;
             int family = 2;
             double dlogc_dnu;
             double *param = new double[2];
             param[0] = rho_double;
             param[1] = nu_double;

             difflPDF_nu_tCopula_new( &u_double, &v_double, &nn, param, &family, &dlogc_dnu);
             delete[] param;
             ops_partials.edge4_.partials_[n] += dlogc_dnu;
         }

      }
      return ops_partials.build(logp);
    }

    template <typename T_u, typename T_v, typename T_rho, typename T_nu>
    inline
    typename stan::return_type<T_u, T_v, T_rho, T_nu>::type
    bicop_student_log(const T_u& u, const T_v& v, const T_rho& rho, const T_nu& nu) {
      return bicop_student_log<false>(u, v, rho, nu);
    }

}
#endif // VIFCOPULA_DISTRIBUTION_BICOP_STUDENT_LOG_HPP


