#ifndef VIFCOPULA_LOGBIFCOP_CPP
#define VIFCOPULA_LOGBIFCOP_CPP

#include <Rcpp.h>
#include <RcppEigen.h>
#include <stan/math.hpp>
#include <omp.h>
#include <bicopdist.h>


// [[Rcpp::depends(RcppEigen)]]
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
        double logLLbicop = 0;
        switch ( copula_type(i,k) ) {
        case 0:
            // Independence copula
                logLLbicop = bicop_independence_log(u,v);
                break;
        case 1:
            // Gaussian copula
                // logLLbicop = vifcopula::bicop_normal_log(u,v,par);
                break;

        default:
                // Code to execute if <variable> does not equal the value following any of the cases
                    break;
        }
        Rcpp::Rcout << " I am here" << std::endl;
        return logLLbicop;
    }
}
#endif // VIFCOPULA_LOGBIFCOP_CPP
