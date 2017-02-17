#ifndef VIFCOPULA_LOGBIFCOP_CPP
#define VIFCOPULA_LOGBIFCOP_CPP

#include <Rcpp.h>
#include <stan/math.hpp>
#include <omp.h>
#include <dist/bicop_independence_log.cpp>
#include <dist/bicop_normal_log.cpp>
#include <dist/bicop_student_log.cpp>
#include <dist/bicop_clayton_log.cpp>
#include <dist/bicop_gumbel_log.cpp>
#include <dist/bicop_frank_log.cpp>
#include <dist/bicop_joe_log.cpp>

// [[Rcpp::depends(StanHeaders)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(openmp)]]

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

    using namespace Rcpp;
    using namespace Eigen;
    using namespace stan::math;
    using namespace vifcopula;

    typedef Eigen::Matrix<int,Eigen::Dynamic,1> vector_int;
    typedef Eigen::Matrix<int,1,Eigen::Dynamic> row_vector_int;
    typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> matrix_int;
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
    typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;

    double logBifcop(const matrix_d& u, const matrix_d& v, const matrix_d& par,
                    const matrix_int& copula_type, const int& t_max,const int& i,
                    const int& k) {
        static const char* function = "vifcopula::logBifcop";
        stan::math::var logLLbicop = 0;
        vector_d u_temp = u.col(i);
        switch ( copula_type(i,k-1) ) { // Note important k.
        case 0:
            // Independence copula
                logLLbicop = bicop_independence_log(u_temp,v);
                Rcpp::Rcout << copula_type(i,k-1) <<  " I am here " << logLLbicop.val() << std::endl;
                break;
        case 1:
            // Gaussian copula
                logLLbicop = vifcopula::bicop_normal_log(u_temp,v,0.5);
                Rcpp::Rcout << copula_type(i,k-1) << " I am here " << logLLbicop.val()  << std::endl;
                break;
        case 2:
            // Student copula
            logLLbicop = vifcopula::bicop_student_log(u_temp,v,0.5, 5);
            Rcpp::Rcout << copula_type(i,k-1) << " I am here " << logLLbicop.val()  << std::endl;
            break;
        case 3:
            // Clayon copula
            logLLbicop = vifcopula::bicop_clayton_log(u_temp,v,0.5);
            Rcpp::Rcout << copula_type(i,k-1) << " I am here " << logLLbicop.val()  << std::endl;
            break;
        case 4:
            // Gumbel copula
            logLLbicop = vifcopula::bicop_gumbel_log(u_temp,v,2);
            Rcpp::Rcout << copula_type(i,k-1) << " I am here " << logLLbicop.val()  << std::endl;
            break;
        case 5:
            // Frank copula
            logLLbicop = vifcopula::bicop_frank_log(u_temp,v,2);
            Rcpp::Rcout << copula_type(i,k-1) << " I am here " << logLLbicop.val()  << std::endl;
            break;
        case 6:
            // Joe copula
            logLLbicop = vifcopula::bicop_joe_log(u_temp,v,3);
            Rcpp::Rcout << copula_type(i,k-1) << " I am here " << logLLbicop.val()  << std::endl;
            break;
        default:
                // Code to execute if <variable> does not equal the value following any of the cases
                    break;
        }
        return logLLbicop.val();
    }
}
#endif // VIFCOPULA_LOGBIFCOP_CPP
