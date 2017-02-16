#ifndef VIFCOPULA_FACTOR_MODEL_CPP
#define VIFCOPULA_FACTOR_MODEL_CPP

#include <Rcpp.h>
#include <stan/math.hpp>
#include <omp.h>
#include <logBifcop.cpp>



// [[Rcpp::depends(StanHeaders)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(openmp)]]


namespace vifcopula {

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


class factor_model{
private:
    int t_max;
    int n_max;
    int k;
    matrix_d u;

    matrix_d v;
    matrix_int copula_type;
    matrix_d par;

    vector_d log_bifcop;

public:
    factor_model(const matrix_d& u_, const matrix_d& v_, const matrix_d& par_,
                 const matrix_int& copula_type_, const int& t_max_,const int& n_max_,
                 const int& k_)
        : u(u_), v(v_),par(par_),
          copula_type(copula_type_), t_max(t_max_), n_max(n_max_),
          k(k_) {
        static const char* function = "vifcopula::factor_model";
        log_bifcop.setZero(n_max);
    }

    void set_v (const matrix_d& v_) {
        static const char* function = "vifcopula::factor_model::set_v";
        v = v_;
    }

    void set_par (const matrix_d& par_) {
        static const char* function = "vifcopula::factor_model::set_rho";
        par = par_;
    }

    void set_copula_type (const matrix_int& copula_type_) {
        static const char* function = "vifcopula::factor_model::set_rho";
        copula_type = copula_type_;
    }

    double log_prob (void) {
        static const char* function = "vifcopula::factor_model::log_prob";

        for (int i = 0; i < n_max; i++){
            log_bifcop[i] = logBifcop(u,v,par,copula_type,t_max,i,k);
        }
        return (log_bifcop.sum());
    }


}; // model

} // namespace
#endif // VIFCOPULA_FACTOR_MODEL_CPP
