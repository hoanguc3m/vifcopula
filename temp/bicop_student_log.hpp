#ifndef VIFCOPULA_DISTRIBUTION_BICOP_STUDENT_LOG_HPP
#define VIFCOPULA_DISTRIBUTION_BICOP_STUDENT_LOG_HPP

#include <stan/math.hpp>
#include <iostream>

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
            const T_partials_return     A_arg = nud2[n];
            const T_partials_return     B_arg = 0.5;

            const T_partials_return     digammaA = digamma(A_arg);
            const T_partials_return     digammaB = digamma(B_arg);
            const T_partials_return     digammaSum = digamma(nud2ph[n]);
            const T_partials_return     betaAB = exp(lbeta(A_arg,B_arg));
            const T_partials_return     nu_p_x1_sq = nu_value[n] + square(inv_u_dbl[n]);
            const T_partials_return     nu_p_x2_sq = nu_value[n] + square(inv_v_dbl[n]);

            const T_partials_return     x_max_int1 = nu_value[n]/ nu_p_x1_sq;
            const T_partials_return     x_max_int2 = nu_value[n]/ nu_p_x2_sq;
            T_partials_return     grad_ibeta_1 = 1;
            T_partials_return     grad_ibeta_2 = 1;

            grad_reg_inc_beta(grad_ibeta_1, grad_ibeta_2, A_arg, B_arg, x_max_int1, digammaA, digammaB, digammaSum, betaAB);
            const T_partials_return dev_ibeta_x1 = grad_ibeta_1;
            grad_reg_inc_beta(grad_ibeta_1, grad_ibeta_2, A_arg, B_arg, x_max_int2, digammaA, digammaB, digammaSum, betaAB);
            const T_partials_return dev_ibeta_x2 = grad_ibeta_1;

            const T_partials_return dx1_dnu = sign(inv_u_dbl[n]) /(2 * pdf(s,inv_u_dbl[n])) * (0.5* dev_ibeta_x1 +
                                                                        pow( nu_p_x1_sq , - nud2ph[n] ) * sign(inv_u_dbl[n]) *
                                                                        pow( nu_value[n], nud2m1[n] ) *  inv_u_dbl[n] / betaAB );

//            std::cout << pdf(s,inv_u_dbl[n]) << " " << A_arg << std::endl;
//            std::cout << B_arg << " " << nud2ph[n] << std::endl;
//            std::cout << pow( nu_value[n], nud2m1[n] ) *  inv_u_dbl[n] << " " <<  pow( nu_p_x1_sq , - nud2ph[n] ) << std::endl;
//            std::cout << betaAB << " " <<  dx1_dnu << std::endl;
            // std::cout << "inbeder_out" << dev_ibeta_x1 << " " << std::endl;
            // std::cout << "inbeder_out" << dev_ibeta_x2 << " " << std::endl;

//              One problem is the imbeder, grad_reg_inc_beta

            const T_partials_return dx2_dnu = sign(inv_v_dbl[n]) /(2 * pdf(s,inv_v_dbl[n])) * (0.5* dev_ibeta_x2 +
                                                                        pow( nu_p_x2_sq , - nud2ph[n] ) * sign(inv_v_dbl[n]) *
                                                                        pow( nu_value[n], nud2m1[n] ) *  inv_v_dbl[n] / betaAB );

            const T_partials_return x1_x1_dnu = 2 * inv_u_dbl[n] * dx1_dnu;
            const T_partials_return x1_x2_dnu = 2 * inv_u_dbl[n] * dx2_dnu;
            const T_partials_return x2_x2_dnu = 2 * inv_v_dbl[n] * dx2_dnu;
            const T_partials_return x2_x1_dnu = 2 * inv_v_dbl[n] * dx1_dnu;


            ops_partials.edge4_.partials_[n] += (- digammaSum) + digammaA + 0.5 * log_1mrhosq[n] -
                                            nud2m1[n] / nu_value[n] - 0.5 * log(nu_value[n]) +
                                            nud2ph[n] * ( (1 + x1_x1_dnu)/nu_p_x1_sq + (1 + x2_x2_dnu)/nu_p_x2_sq ) +
                                            0.5 * ( log(nu_p_x1_sq) + log(nu_p_x2_sq)) -
                                            nud2p1[n] *(1-sq_rho[n] + x1_x1_dnu + x2_x2_dnu - rho_value[n]*(x1_x2_dnu+x2_x1_dnu) )/ M_nu_rho -
                                            0.5 * logM;
//            std::cout << dev_ibeta_x1 << " " << dev_ibeta_x2 << std::endl;
//            std::cout << dx1_dnu << " " << dx2_dnu << std::endl;
//            std::cout << digammaSum << " " << digammaA << std::endl;
//            std::cout << 0.5 * log_1mrhosq[n]  << " " << nud2m1[n] / nu_value[n] << std::endl;
//            std::cout << 0.5 * log(nu_value[n]) << " " << (- digammaSum) + digammaA + 0.5 * log_1mrhosq[n] -
//                                            nud2m1[n] / nu_value[n] - 0.5 * log(nu_value[n]) << std::endl;
//            std::cout << (1 + x1_x1_dnu) << " " << (1 + x2_x2_dnu) << std::endl;
//            std::cout << nud2ph[n] * ( (1 + x1_x1_dnu)/nu_p_x1_sq + (1 + x2_x2_dnu)/nu_p_x2_sq ) << " " << M_nu_rho << std::endl;
//            std::cout << nud2p1[n] << " " << (1-sq_rho[n] + (1 - rho_value[n]) * (x1_x1_dnu + x2_x2_dnu) ) << std::endl;
//            std::cout << 0.5 * ( log(nu_p_x1_sq) + log(nu_p_x2_sq)) << " " << 0.5 * logM << std::endl;

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


