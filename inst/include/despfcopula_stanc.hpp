// Code generated by Stan version 2.14

#ifndef VIFCOPULA_DESP_FACTOR_COP_HPP
#define VIFCOPULA_DESP_FACTOR_COP_HPP

#include <omp.h>
#include <stan/model/model_header.hpp>
#include <dist/bicop_log.hpp>
#include <service/write_theta.hpp>

// [[Rcpp::depends(StanHeaders)]]
// [[Rcpp::depends(rstan)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(openmp)]]


namespace vifcopula
{

using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
using namespace vifcopula;

typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;
typedef Eigen::Matrix<int,Eigen::Dynamic,1> vector_int;
typedef Eigen::Matrix<int,1,Eigen::Dynamic> row_vector_int;
typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> matrix_int;


//static int current_statement_begin__;

class despfcopula : public prob_grad
{
private:
    int t_max;                  // Number of periods
    int n_max;                  // Number of series
    int k;                      // Number of Latents
    matrix_d u;                 // Input matrix copula[t*n]
    matrix_d u_eps;                 // Input matrix copula[t*n]
    vector<int> gid;             // Group of copula
    vector<int> copula_type;     // copula type
    vector<int> copula_eps_type;     // latent copula type
    int twofcop;

public:

    ~despfcopula() { }

    despfcopula(const matrix_d& u_,
                const matrix_d& u_eps_,
                const vector<int>& gid_,
              const vector<int>& copula_type_,
              const vector<int>& copula_eps_type_,
              int& twofcop_,
              const int& t_max_,const int& n_max_,const int& k_,
              std::ostream* pstream__ = 0)
        : prob_grad(0)
    {
        typedef boost::ecuyer1988 rng_t;
        rng_t base_rng(0);  // 0 seed default
        despfcopula(u_, u_eps_,gid_,
            copula_type_, copula_eps_type_, twofcop_,
                  t_max_, n_max_, k_,
                  base_rng, pstream__);
    }

    void recheck_num_para_r(size_t& num_params_r__){

        num_params_r__ = 0U;
        param_ranges_i__.clear();
        num_params_r__ += t_max*2;

        num_params_r__ += 1;
        if (twofcop == 0) {
            num_params_r__ --;
        }
        else if (twofcop == 2) {
            num_params_r__ ++;
        }

        num_params_r__ += n_max;
        for (int i = 0; i < n_max; i++) {
            if (copula_type[i] == 0) {
                num_params_r__ --;
            }
            else if (copula_type[i] == 2) {
                num_params_r__ ++;
            }
        }

        num_params_r__ += n_max;
        for (int i = 0; i < n_max; i++) {
            if (copula_eps_type[i] == 0) {
                num_params_r__ --;
            }
            else if (copula_eps_type[i] == 2) {
                num_params_r__ ++;
            }
        }


    }

    template <class RNG>
    despfcopula(const matrix_d& u_, const matrix_d& u_eps_,
                const vector<int>& gid_,
              const vector<int>& copula_type_,
              const vector<int>& copula_eps_type_,
              const int& twofcop_,
              const int& t_max_,const int& n_max_,const int& k_,
              RNG& base_rng__,std::ostream* pstream__ = 0)
        : prob_grad(0), u(u_), u_eps(u_eps_),gid(gid_),
            copula_type(copula_type_),copula_eps_type(copula_eps_type_),twofcop(twofcop_),
          t_max(t_max_), n_max(n_max_),k(k_)
    {
        static const char* function = "vifcopula::despfcopula";

        recheck_num_para_r(num_params_r__);


    }

    void set_copula_type(std::vector<int> copula_type_)
    {
        copula_type = copula_type_;
        recheck_num_para_r(num_params_r__);
    }

    void set_copula_eps_type(std::vector<int> copula_eps_type_)
    {
        copula_eps_type = copula_eps_type_;
        recheck_num_para_r(num_params_r__);
    }
    void set_copula_twofcop(int twofcop_)
    {
        twofcop = twofcop_;
        recheck_num_para_r(num_params_r__);
    }

    int get_eff_para(void){
        int eff_para = num_params_r__ - t_max*2;
        return eff_para;
    }
    int get_t_max(void){
        return t_max;
    }
    double calc_log_over_v( RNG& base_rng__,
                            Eigen::VectorXd& mean_vi,
                            int eff_num_para){

        return 0;
    }

    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const
    {

        T__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        static int current_statement_begin__;

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        // model parameters
        stan::io::reader<T__> in__(params_r__,params_i__);

        Eigen::Matrix<T__,Eigen::Dynamic,1>  v_d;
        (void) v_d;  // dummy to suppress unused var warning
        if (jacobian__)
            v_d = in__.vector_lub_constrain(0,1,t_max,lp__);
        else
            v_d = in__.vector_lub_constrain(0,1,t_max);

        Eigen::Matrix<T__,Eigen::Dynamic,1>  v_eps;
        (void) v_eps;  // dummy to suppress unused var warning
        if (jacobian__)
            v_eps = in__.vector_lub_constrain(0,1,t_max, lp__);
        else
            v_eps = in__.vector_lub_constrain(0,1,t_max);

        Eigen::Matrix<T__,Eigen::Dynamic,1>  theta_twofcop(1);
        (void) theta_twofcop;  // dummy to suppress unused var warning
        Eigen::Matrix<T__,Eigen::Dynamic,1>  theta2_twofcop(1);
        (void) theta2_twofcop;  // dummy to suppress unused var warning

        Eigen::Matrix<T__,Eigen::Dynamic,1>  theta(n_max);
        (void) theta;  // dummy to suppress unused var warning

        Eigen::Matrix<T__,Eigen::Dynamic,1>  theta2(n_max);
        (void) theta2;  // dummy to suppress unused var warning


        Eigen::Matrix<T__,Eigen::Dynamic,1>  theta_eps(n_max);
        (void) theta_eps;  // dummy to suppress unused var warning

        Eigen::Matrix<T__,Eigen::Dynamic,1>  theta2_eps(n_max);
        (void) theta2_eps;  // dummy to suppress unused var warning


        vector_d u_col;
        // Eigen::Matrix<T__,Eigen::Dynamic,1> v_d_col;
        // Eigen::Matrix<T__,Eigen::Dynamic,1> v_eps_col;

        // model body
        try
        {

            current_statement_begin__ = 12;

            std::vector<int> twofcop_vec(1);
            twofcop_vec[0] = twofcop;
            bicop_log_add<propto__,jacobian__,T__,T__,T__>(0, twofcop_vec, v_d, v_eps, theta_twofcop, theta2_twofcop, lp__, lp_accum__, in__);


            for (int i = 0; i < n_max; i++) {
                u_col = u.col(i);
                bicop_log_add<propto__,jacobian__,double,T__,T__>(i, copula_type, u_col, v_d, theta, theta2, lp__, lp_accum__, in__);
            }

            current_statement_begin__ = 13;

            for (int i = 0; i < n_max; i++) {
                u_col = u_eps.col(i);
                bicop_log_add<propto__,jacobian__,double,T__,T__>(i, copula_eps_type, u_col, v_eps, theta_eps, theta2_eps, lp__, lp_accum__, in__);
            }
        }
        catch (const std::exception& e)
        {
            stan::lang::rethrow_located(e,current_statement_begin__);
            throw std::runtime_error(" ");
        }

        lp_accum__.add(lp__);
        return lp_accum__.sum();

    } // log_prob()

    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                std::ostream* pstream = 0) const
    {
        std::vector<T_> vec_params_r;
        vec_params_r.reserve(params_r.size());
        for (int i = 0; i < params_r.size(); ++i)
            vec_params_r.push_back(params_r(i));
        std::vector<int> vec_params_i;
        return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }


    void get_param_names(std::vector<std::string>& names__) const
    {
        names__.resize(0);
        names__.push_back("v");
        names__.push_back("theta");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const
    {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(t_max*2);
        dims__.push_back(num_params_r__ - t_max*2);
    }

    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const
    {
        vars__.resize(0);
        stan::io::reader<double> in__(params_r__,params_i__);
        static const char* function__ = "vifcopula::write_array";
        (void) function__; // dummy call to supress warning
        // read-transform, write parameters
        vector_d v = in__.vector_lub_constrain(0,1,t_max*2);

        for (int k_0__ = 0; k_0__ < t_max*2; ++k_0__)
        {
            vars__.push_back(v[k_0__]);
        }

        int num_theta_param = num_params_r__ -  t_max*2;

        write_theta(0, twofcop, in__, vars__);

        for (int i = 0; i < n_max; i++) {
            write_theta(i, copula_type[i], in__, vars__);
        }

        for (int i = 0; i < n_max; i++) {
            write_theta(i, copula_eps_type[i], in__, vars__);
        }

        if (!include_tparams__) return;
        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__; // dummy call to supress warning
        stan::math::accumulator<double> lp_accum__;

        double DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning



        // validate transformed parameters

        // write transformed parameters

        if (!include_gqs__) return;
        // declare and define generated quantities

    }

    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const
    {
        std::vector<double> params_r_vec(params_r.size());
        for (int i = 0; i < params_r.size(); ++i)
            params_r_vec[i] = params_r(i);
        std::vector<double> vars_vec;
        std::vector<int> params_i_vec;
        write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
        vars.resize(vars_vec.size());
        for (int i = 0; i < vars.size(); ++i)
            vars(i) = vars_vec[i];
    }

    static std::string model_name()
    {
        return "despfcopula";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const
    {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= t_max*2; ++k_0__)
        {
            param_name_stream__.str(std::string());
            param_name_stream__ << "v" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= (num_params_r__ - t_max*2); ++k_0__)
        {
            param_name_stream__.str(std::string());
            param_name_stream__ << "theta" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const
    {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= t_max*2; ++k_0__)
        {
            param_name_stream__.str(std::string());
            param_name_stream__ << "v" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= (num_params_r__ - t_max*2); ++k_0__)
        {
            param_name_stream__.str(std::string());
            param_name_stream__ << "theta" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
    }

}; // model

} // namespace


#endif // VIFCOPULA_DESP_FACTOR_COP_HPP

