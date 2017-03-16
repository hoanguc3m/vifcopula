#ifndef VIFCOPULA_ONE_FACTOR_COP_CPP
#define VIFCOPULA_ONE_FACTOR_COP_CPP

#include <Rcpp.h>
//#include <rstan/rstaninc.hpp>
#include <stan/model/model_header.hpp>
#include <stan/math.hpp>

#include <omp.h>
#include <logBifcop.cpp>



// [[Rcpp::depends(StanHeaders)]]
// [[Rcpp::depends(rstan)]]
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


class one_factor_cop : public prob_grad {
private:
    int t_max;
    int n_max;
    int k;
    matrix_d u;
    vector_int gid;
    vector_int copula_type;

public:
    one_factor_cop(const matrix_d& u_, const vector_int& gid_,
                   const matrix_int& copula_type_
                   const int& t_max_,const int& n_max_,const int& k_)
        : u(u_), gid(gid_),copula_type(copula_type_),
          t_max(t_max_), n_max(n_max_),k(k_) {
        static const char* function = "vifcopula::one_factor_cop";

        // set parameter ranges
        num_params_r__ = 0U;
        param_ranges_i__.clear();
        num_params_r__ += t_max;
        num_params_r__ += n_max;

    }

    // void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
    //     // get dim based on copula type
    //     dimss__.resize(0);
    //     std::vector<size_t> dims__;
    //     dims__.resize(0);
    //     dims__.push_back(t_max);
    //     dimss__.push_back(dims__);
    //
    //     dims__.resize(0);
    //     dims__.push_back(n_max);
    //     dimss__.push_back(dims__);
    // }

    void set_copula_type (const vector_int& copula_type_) {
        static const char* function = "vifcopula::one_factor_cop::set_copula_type";
        copula_type = copula_type_;
    }

    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        static const char* function = "vifcopula::one_factor_cop::log_prob";

        vector_d log_bifcop;
        log_bifcop.setZero(n_max);

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        // model parameters
        stan::io::reader<T__> in__(params_r__,params_i__);

        Eigen::Matrix<T__,Eigen::Dynamic,1>  v;
        (void) v;  // dummy to suppress unused var warning
        if (jacobian__)
            v = in__.vector_lub_constrain(-(1),1,t_max,lp__);
        else
            v = in__.vector_lub_constrain(-(1),1,t_max);

        Eigen::Matrix<T__,Eigen::Dynamic,1>  par;
        (void) par;  // dummy to suppress unused var warning
        if (jacobian__)
            par = in__.vector_lub_constrain(-(1),1,n_max,lp__);
        else
            par = in__.vector_lub_constrain(-(1),1,n_max);

        // model body
        try {

            current_statement_begin__ = 12;
            lp_accum__.add(uniform_log<propto__>(v, -(1), 1));
            current_statement_begin__ = 13;
            lp_accum__.add(uniform_log<propto__>(par, -(1), 1));
            current_statement_begin__ = 14;

            int k_temp = k-1;
            for (int i = 0; i < n_max; i++){
                current_statement_begin__ = 15;
                lp_accum__.add(logBifcop(u,v,par,copula_type,t_max,i,k_temp));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e,current_statement_begin__);
            // Next line prevents compiler griping about no return
            throw std::runtime_error(" ");
        }

        lp_accum__.add(lp__);
        return lp_accum__.sum();
    }

    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                std::ostream* pstream = 0) const {
        std::vector<T_> vec_params_r;
        vec_params_r.reserve(params_r.size());
        for (int i = 0; i < params_r.size(); ++i)
            vec_params_r.push_back(params_r(i));
        std::vector<int> vec_params_i;
        return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }

    ~ one_factor_cop(){}


}; // model

} // namespace
#endif // VIFCOPULA_ONE_FACTOR_COP_CPP
