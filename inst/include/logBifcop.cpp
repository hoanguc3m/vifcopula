#ifndef VIFCOPULA_LOGBIFCOP_CPP
#define VIFCOPULA_LOGBIFCOP_CPP

#include <stan/math.hpp>
#include <dist/bicop_independence_log.hpp>
#include <dist/bicop_normal_log.hpp>
#include <dist/bicop_student_log.hpp>
#include <dist/bicop_clayton_log.hpp>
#include <dist/bicop_gumbel_log.hpp>
#include <dist/bicop_frank_log.hpp>
#include <dist/bicop_joe_log.hpp>

namespace vifcopula {
/**
 * The log of the bivariate copula density for the specified type copula_type, specified vector(s) u
 * and the specified vector(s) v given parameter(s) theta, theta2. u, v, or theta, theta2 can
 * each be either a scalar or a vector. Any vector inputs
 * must be the same length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation/observation/theta/theta2.
 * @param u (Sequence of) scalar(s) in [0,1].
 * @param v (Sequence of) scalar(s) in [0,1].
 * @param theta (Sequence of) first parameters for the bivariate copula
 * @param theta2 (Sequence of) second parameters for the bivariate copula
 * @return The log of the product of the densities.
 * @throw std::domain_error if the scale is not positive.
 * @tparam T_u Underlying type of scalar in sequence.
 * @tparam T_v Underlying type of scalar in sequence.
 * @tparam T_theta Type of first parameter.
 * @tparam T_theta2 Type of second parameter.
 */

    using namespace stan::math;
    using namespace vifcopula;

    template <bool propto, typename T_u, typename T_v, typename T_theta, typename T_theta2>
    typename stan::return_type<T_u, T_v, T_theta, T_theta2>::type
        logBifcop(const int& copula_type, const T_u& u, const T_v& v, const T_theta& theta = 0, const T_theta2& theta2 = 0) {
            static const char* function("vifcopula::logBifcop");

            var logLLbicop; // Note find the return type.

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


    template <  typename T_u,
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
