// Code generated by Stan version 2.14

#ifndef VIFCOPULA_ONE_FACTOR_COP_HPP
#define VIFCOPULA_ONE_FACTOR_COP_HPP

#include <omp.h>
#include <stan/model/model_header.hpp>
#include <dist/bicop_log.hpp>
#include <service/write_theta.hpp>

#include <boost/math/quadrature/gauss_kronrod.hpp>


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

class onefcopula : public prob_grad
{
private:
    int t_max;                  // Number of periods
    int n_max;                  // Number of series
    int k;                      // Latent index
    matrix_d u;                 // Input matrix copula[t*n]
    vector<int> gid;             // Group of copula
    vector<int> copula_type;     // copula type

public:
    /* Generated by mcstan 2.14
    onefcopula(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        typedef boost::ecuyer1988 rng_t;
        rng_t base_rng(0);  // 0 seed default
        ctor_body(context__, base_rng, pstream__);
    }
    template <class RNG>
    onefcopula(stan::io::var_context& context__,
        RNG& base_rng__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, base_rng__, pstream__);
    }

    template <class RNG>
    void ctor_body(stan::io::var_context& context__,
                   RNG& base_rng__,
                   std::ostream* pstream__) {
        current_statement_begin__ = -1;

        static const char* function__ = "vifcopula::onefcopula";
        (void) function__; // dummy call to supress warning
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        double DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        // initialize member variables
        context__.validate_dims("data initialization", "n_max", "int", context__.to_vec());
        n_max = int(0);
        vals_i__ = context__.vals_i("n_max");
        pos__ = 0;
        n_max = vals_i__[pos__++];
        context__.validate_dims("data initialization", "t_max", "int", context__.to_vec());
        t_max = int(0);
        vals_i__ = context__.vals_i("t_max");
        pos__ = 0;
        t_max = vals_i__[pos__++];
        context__.validate_dims("data initialization", "u", "matrix_d", context__.to_vec(t_max,n_max));
        validate_non_negative_index("u", "t_max", t_max);
        validate_non_negative_index("u", "n_max", n_max);
        u = matrix_d(static_cast<Eigen::VectorXd::Index>(t_max),static_cast<Eigen::VectorXd::Index>(n_max));
        vals_r__ = context__.vals_r("u");
        pos__ = 0;
        size_t u_m_mat_lim__ = t_max;
        size_t u_n_mat_lim__ = n_max;
        for (size_t n_mat__ = 0; n_mat__ < u_n_mat_lim__; ++n_mat__) {
            for (size_t m_mat__ = 0; m_mat__ < u_m_mat_lim__; ++m_mat__) {
                u(m_mat__,n_mat__) = vals_r__[pos__++];
            }
        }
        context__.validate_dims("data initialization", "k", "int", context__.to_vec());
        k = int(0);
        vals_i__ = context__.vals_i("k");
        pos__ = 0;
        k = vals_i__[pos__++];
        context__.validate_dims("data initialization", "gid", "int", context__.to_vec(n_max));
        validate_non_negative_index("gid", "n_max", n_max);
        gid = std::vector<int>(n_max,int(0));
        vals_i__ = context__.vals_i("gid");
        pos__ = 0;
        size_t gid_limit_0__ = n_max;
        for (size_t i_0__ = 0; i_0__ < gid_limit_0__; ++i_0__) {
            gid[i_0__] = vals_i__[pos__++];
        }
        context__.validate_dims("data initialization", "copula_type", "int", context__.to_vec(n_max));
        validate_non_negative_index("copula_type", "n_max", n_max);
        copula_type = std::vector<int>(n_max,int(0));
        vals_i__ = context__.vals_i("copula_type");
        pos__ = 0;
        size_t copula_type_limit_0__ = n_max;
        for (size_t i_0__ = 0; i_0__ < copula_type_limit_0__; ++i_0__) {
            copula_type[i_0__] = vals_i__[pos__++];
        }

        // validate, data variables
        check_greater_or_equal(function__,"n_max",n_max,0);
        check_greater_or_equal(function__,"t_max",t_max,0);
        check_greater_or_equal(function__,"k",k,0);
        // initialize data variables

        try {
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e,current_statement_begin__);
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }

        // validate transformed data

        // set parameter ranges
        num_params_r__ = 0U;
        param_ranges_i__.clear();
        num_params_r__ += t_max;
        num_params_r__ += n_max;
        for (int i = 0; i < n_max; i++){
            if (copula_type[i] == 0){
                num_params_r__ --;
            } else if (copula_type[i] == 2){
                num_params_r__ ++;
            }

        }
    }

        void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;

        if (!(context__.contains_r("v")))
            throw std::runtime_error("variable v missing");
        vals_r__ = context__.vals_r("v");
        pos__ = 0U;
        context__.validate_dims("initialization", "v", "vector_d", context__.to_vec(t_max));
        // generate_declaration v
        vector_d v(static_cast<Eigen::VectorXd::Index>(t_max));
        for (int j1__ = 0U; j1__ < t_max; ++j1__)
            v(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_lub_unconstrain(0,1,v);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Error transforming variable v: ") + e.what());
        }

        if (!(context__.contains_r("theta")))
            throw std::runtime_error("variable theta missing");
        vals_r__ = context__.vals_r("theta");
        pos__ = 0U;
        context__.validate_dims("initialization", "theta", "vector_d", context__.to_vec(n_max));
        // generate_declaration theta
        vector_d theta(static_cast<Eigen::VectorXd::Index>(n_max));
        for (int j1__ = 0U; j1__ < n_max; ++j1__)
            theta(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_lub_unconstrain(-(1),1,theta);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Error transforming variable theta: ") + e.what());
        }

        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    */

    ~onefcopula() { }

    onefcopula(const matrix_d& u_, const vector<int>& gid_,
               const vector<int>& copula_type_,
               const int& t_max_,const int& n_max_,const int& k_,
               std::ostream* pstream__ = 0)
        : prob_grad(0)
    {
        typedef boost::ecuyer1988 rng_t;
        rng_t base_rng(0);  // 0 seed default
        onefcopula(u_, gid_, copula_type_,
                   t_max_, n_max_, k_,
                   base_rng, pstream__);
    }

    template <class RNG>
    onefcopula(const matrix_d& u_, const vector<int>& gid_,
               const vector<int>& copula_type_,
               const int& t_max_,const int& n_max_,const int& k_,
               RNG& base_rng__,std::ostream* pstream__ = 0)
        : prob_grad(0), u(u_), gid(gid_),copula_type(copula_type_),
          t_max(t_max_), n_max(n_max_),k(k_)
    {
        static const char* function = "vifcopula::onefcopula";

        // Inherit from prob_grad
        //      size_t num_params_r__;
        //      std::vector<std::pair<int, int> > param_ranges_i__;
        // set parameter ranges
        num_params_r__ = 0U;
        param_ranges_i__.clear();
        num_params_r__ += t_max;
        num_params_r__ += n_max;
        for (int i = 0; i < n_max; i++) {
            if (copula_type[i] == 0) {
                num_params_r__ --;
            }
            else if (copula_type[i] == 2) {
                num_params_r__ ++;
            }
        }


    }

    void set_copula_type(std::vector<int> copula_type_) {
        copula_type = copula_type_;

        num_params_r__ = 0U;
        param_ranges_i__.clear();
        num_params_r__ += t_max;
        num_params_r__ += n_max;
        for (int i = 0; i < n_max; i++) {
            if (copula_type[i] == 0) {
                num_params_r__ --;
            }
            else if (copula_type[i] == 2) {
                num_params_r__ ++;
            }
        }
    }

    int get_eff_para(void){
        int eff_para = 0;
        eff_para += n_max;
        for (int i = 0; i < n_max; i++) {
            if (copula_type[i] == 0) {
                eff_para --;
            }
            else if (copula_type[i] == 2) {
                eff_para ++;
            }
        }
        return eff_para;
    }
    int get_t_max(void){
        return t_max;
    }
    int get_n_max(void){
        return n_max;
    }
    // template <typename RNG>
    // double calc_log_over_v( RNG& base_rng__,
    //                         Eigen::VectorXd& mean_vi,
    //                         int eff_num_para){
    //     std::srand(base_rng__());
    //
    //     double logc = 0;
    //     int MCnum = 1000;
    //     vector<double> logc_t(t_max,0.0);
    //
    //     Eigen::VectorXd theta_12 = mean_vi.tail(eff_num_para);
    //     int count = 0;
    //     vector<double> theta(n_max,0.0);
    //     vector<double> theta2(n_max,0.0);
    //
    //     for (int i = 0; i < n_max; i++) {
    //         if (copula_type[i] == 2) {
    //             theta[i] = theta_12(count); count++;
    //             theta2[i] = theta_12(count); count++;
    //         } else {
    //             theta[i] = theta_12(count); count++;
    //         }
    //     }
    //     for (int t = 0; t < t_max; t++) {
    //         Eigen::VectorXd v_t = (Eigen::VectorXd::Random(MCnum)).array().abs() ; // range [-1,1]
    //         Eigen::VectorXd logc_jt = Eigen::VectorXd::Zero(MCnum);
    //
    //         for (int j = 0; j < MCnum; j++) {
    //             for (int i = 0; i < n_max; i++) {
    //                 logc_jt(j) += bicop_log_double(copula_type[i], u(t,i), v_t(j), theta[i], theta2[i] )   ;
    //             }
    //         }
    //
    //         double max_logct = logc_jt.maxCoeff();
    //         logc_t[t] = max_logct + log ( (logc_jt.array() - max_logct).array().exp().sum())   ;
    //         logc_t[t] -= log(MCnum);
    //     }
    //
    //     for (auto& log_val : logc_t)
    //         logc += log_val;
    //     return logc;
    // }

    template <typename RNG>
    double calc_log_over_v( RNG& base_rng__,
                            Eigen::VectorXd& mean_vi,
                            int eff_num_para){
        std::srand(base_rng__());

        double logc = 0;
        const int MCnum = 25;
        const ulong MCnum_n = MCnum / 2 + 1;

        const std::array<double, MCnum_n> ab = boost::math::quadrature::gauss<double, MCnum>::abscissa();
        const std::array<double, MCnum_n> w = boost::math::quadrature::gauss<double, MCnum>::weights();

        Eigen::VectorXd v_t(MCnum);
        Eigen::VectorXd weight_t(MCnum);

        v_t[0] = ab[0] *0.5 + 0.5; // Tranform from [-1,1]
        weight_t[0] = w[0] * 0.5;// Tranform from [-1,1] to [0,1]

        for (unsigned i = 1; i < ab.size(); ++i){
            v_t[2*i-1] = ab[i] *0.5 + 0.5;
            v_t[2*i] = - ab[i] *0.5 + 0.5;
            weight_t[2*i-1] = w[i] * 0.5;
            weight_t[2*i] = w[i] * 0.5;
        }


        vector<double> logc_t(t_max,0.0);

        Eigen::VectorXd theta_12 = mean_vi.tail(eff_num_para);
        int count = 0;
        vector<double> theta(n_max,0.0);
        vector<double> theta2(n_max,0.0);

        for (int i = 0; i < n_max; i++) {
            if (copula_type[i] == 2) {
                theta[i] = theta_12(count); count++;
                theta2[i] = theta_12(count); count++;
            } else {
                if (copula_type[i] != 0) {
                    theta[i] = theta_12(count); count++;
                }

            }
        }
        for (int t = 0; t < t_max; t++) {
            Eigen::VectorXd logc_jt = Eigen::VectorXd::Zero(MCnum);

            for (int j = 0; j < MCnum; j++) {
                for (int i = 0; i < n_max; i++) {
                    logc_jt(j) += bicop_log_double(copula_type[i], u(t,i), v_t(j), theta[i], theta2[i] )   ;
                }
            }

            double max_logct = logc_jt.maxCoeff();
            Eigen::VectorXd exp_logc_jt_minus_max = (logc_jt.array() - max_logct).array().exp();
            logc_t[t] = max_logct + log ( (exp_logc_jt_minus_max.array() * weight_t.array() ).sum())   ;

        }

        for (auto& log_val : logc_t){
            logc += log_val;
        }

        return logc;
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

        double UMAX = 1-1e-10;
        double UMIN = 1e-10;
        Eigen::Matrix<T__,Eigen::Dynamic,1>  v;
        (void) v;  // dummy to suppress unused var warning
        if (jacobian__)
            v = in__.vector_lub_constrain(UMIN,UMAX,t_max,lp__);
        else
            v = in__.vector_lub_constrain(UMIN,UMAX,t_max);

        Eigen::Matrix<T__,Eigen::Dynamic,1>  theta(n_max);
        (void) theta;  // dummy to suppress unused var warning

        Eigen::Matrix<T__,Eigen::Dynamic,1>  theta2(n_max);
        (void) theta2;  // dummy to suppress unused var warning

        vector_d u_col;

        // model body
        try {

            current_statement_begin__ = 12;

            // if (include_summand<propto__>::value)
            //     lp_accum__.add(uniform_lpdf<propto__>(v, -(1), 1));

            current_statement_begin__ = 13;


            for (int i = 0; i < n_max; i++) {
                current_statement_begin__ = 14;
                u_col = u.col(i);

                bicop_log_add<propto__,jacobian__,double, T__, T__>(i, copula_type, u_col, v, theta, theta2, lp__, lp_accum__, in__);


            }
        }
        catch (const std::exception& e)
        {
            stan::lang::rethrow_located(e,current_statement_begin__);
            // Next line prevents compiler griping about no return
            throw std::runtime_error(" ");
        }
        //std::cout << " End_of_log_prob " << theta << std::endl;

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
        dims__.push_back(t_max);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(num_params_r__ - t_max);
        dimss__.push_back(dims__);
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
        double UMAX = 1-1e-10;
        double UMIN = 1e-10;
        vector_d v = in__.vector_lub_constrain(UMIN,UMAX,t_max);

        for (int k_0__ = 0; k_0__ < t_max; ++k_0__)
        {
            vars__.push_back(v[k_0__]);
        }

        /* Generated by mcstan 2.14
        vector_d theta = in__.vector_lub_constrain(-(1),1,n_max);
        for (int k_0__ = 0; k_0__ < n_max; ++k_0__) {
            vars__.push_back(theta[k_0__]);
        }
        */
        int num_theta_param = num_params_r__ -  t_max;

        for (int i = 0; i < n_max; i++) {
            write_theta(i,copula_type[i], in__, vars__);
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


        // try {
        // } catch (const std::exception& e) {
        //     stan::lang::rethrow_located(e,current_statement_begin__);
        //     // Next line prevents compiler griping about no return
        //     throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        // }

        // validate generated quantities

        // write generated quantities
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
        return "onefcopula";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const
    {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= t_max; ++k_0__)
        {
            param_name_stream__.str(std::string());
            param_name_stream__ << "v" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= (num_params_r__ - t_max); ++k_0__)
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
        for (int k_0__ = 1; k_0__ <= t_max; ++k_0__)
        {
            param_name_stream__.str(std::string());
            param_name_stream__ << "v" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        for (int k_0__ = 1; k_0__ <= (num_params_r__ - t_max); ++k_0__)
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

//typedef vifcopula::onefcopula stan_model;
#endif // VIFCOPULA_ONE_FACTOR_COP_HPP

