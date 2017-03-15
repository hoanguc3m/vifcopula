#ifndef VIFCOPULA_LOGBIFCOP_CPP
#define VIFCOPULA_LOGBIFCOP_CPP

#include <stan/math.hpp>
#include <dist/bicop_independence_log.cpp>
#include <dist/bicop_normal_log.cpp>
#include <dist/bicop_student_log.cpp>
#include <dist/bicop_clayton_log.cpp>
#include <dist/bicop_gumbel_log.cpp>
#include <dist/bicop_frank_log.cpp>
#include <dist/bicop_joe_log.cpp>

namespace vifcopula {
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
 * @throw std::domain_error if the scale is not positive.
 * @tparam T_u Underlying type of scalar in sequence.
 * @tparam T_v Underlying type of scalar in sequence.
 * @tparam T_rho Type of correlation parameter.
 */

    using namespace stan::math;
    using namespace vifcopula;

    template <bool propto, typename T_coptype, typename T_u, typename T_v, typename T_theta, typename T_theta2>
    typename stan::return_type<T_u, T_v, T_theta, T_theta2>::type
        logBifcop(const int& copula_type, const T_u& u, const T_v& v, const T_theta& theta = 0, const T_theta2& theta2 = 0) {
            static const char* function("vifcopula::logBifcop");

            stan::return_type<T_u, T_v, T_theta, T_theta2> logLLbicop = 0; // Note find the return type.

            switch ( copula_type ) {
            case 0:
                // Independence copula
                logLLbicop = bicop_independence_log(u,v);
                break;
            case 1:
                // Gaussian copula
                logLLbicop = vifcopula::bicop_normal_log(u,v,theta);
                break;
            case 2:
                // Student copula
                logLLbicop = vifcopula::bicop_student_log(u,v,theta, theta2);
                break;
            case 3:
                // Clayon copula
                logLLbicop = vifcopula::bicop_clayton_log(u,v,theta);
                break;
            case 4:
                // Gumbel copula
                logLLbicop = vifcopula::bicop_gumbel_log(u,v,theta);
                break;
            case 5:
                // Frank copula
                logLLbicop = vifcopula::bicop_frank_log(u,v,theta);
                break;
            case 6:
                // Joe copula
                logLLbicop = vifcopula::bicop_joe_log(u,v,theta);
                break;
            default:
                // Code to execute if <variable> does not equal the value following any of the cases
                // Send a break message.
                break;
            }
            //        Rcpp::Rcout << copula_type(i,k) << " I am here " << logLLbicop.val()  << std::endl;
            return logLLbicop;
        }


    template <typename T_coptype,
                typename T_u,
                typename T_v,
                typename T_theta = double,
                typename T_theta2 = double>
    inline
        typename stan::return_type<T_u, T_v, T_theta, T_theta2>::type
        logBifcop(const int& copula_type, const T_u& u, const T_v& v, const T_theta& theta = 0, const T_theta2& theta2 = 0) {
            return logBifcop<false>(copula_type, u, v, theta,theta2);
    }



}
#endif // VIFCOPULA_LOGBIFCOP_CPP
