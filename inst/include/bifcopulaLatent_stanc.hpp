// Code generated by Stan version 2.14

#ifndef VIFCOPULA_BI_FACTOR_COP_LATENT_HPP
#define VIFCOPULA_BI_FACTOR_COP_LATENT_HPP

#include <omp.h>
#include <stan/model/model_header.hpp>
#include <dist/bicop_log.hpp>
#include <transform/hfunc_stan.hpp>
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

class bifcopulaLatent : public prob_grad
{
private:
    int t_max;                  // Number of periods
    int n_max;                  // Number of series
    int k;                      // Number of Latents
    matrix_d u;                 // Input matrix copula[t*n]
    vector_d theta;      //parameter1
    vector_d theta2;     //parameter2
    vector_d latent_theta;      //latent_theta1
    vector_d latent_theta2;     //latent_theta2

    vector<int> gid;             // Group of copula
    vector<int> copula_type;     // copula type
    vector<int> latent_copula_type;     // latent copula type

public:

    ~bifcopulaLatent() { }

    bifcopulaLatent(const matrix_d& u_,
            const vector_d& theta_,
            const vector_d& theta2_,
            const vector_d& latent_theta_,
            const vector_d& latent_theta2_,
              vector<int>& gid_,
              const vector<int>& copula_type_,
              const vector<int>& latent_copula_type_,
              const int& t_max_,const int& n_max_,const int& k_,
              std::ostream* pstream__ = 0)
        : prob_grad(0)
    {
        typedef boost::ecuyer1988 rng_t;
        rng_t base_rng(0);  // 0 seed default
        bifcopulaLatent(u_, theta_, theta2_,
                        latent_theta_, latent_theta2_,
                  gid_, copula_type_,
                  latent_copula_type_,
                  t_max_, n_max_, k_,
                  base_rng, pstream__);
    }

    template <class RNG>
    bifcopulaLatent(const matrix_d& u_,
                    const vector_d& theta_,
                    const vector_d& theta2_,
                    const vector_d& latent_theta_,
                    const vector_d& latent_theta2_,
                    const vector<int>& gid_,
              const vector<int>& copula_type_,
              const vector<int>& latent_copula_type_,
              const int& t_max_,const int& n_max_,const int& k_,
              RNG& base_rng__,std::ostream* pstream__ = 0)
        : prob_grad(0), u(u_), theta(theta_), theta2(theta2_),
          latent_theta(latent_theta_), latent_theta2(latent_theta2_),
          gid(gid_),copula_type(copula_type_),
          latent_copula_type(latent_copula_type_),
          t_max(t_max_), n_max(n_max_),k(k_)
    {
        static const char* function = "vifcopula::bifcopulaLatent";

        // Inherit from prob_grad
        //      size_t num_params_r__;
        //      std::vector<std::pair<int, int> > param_ranges_i__;
        // set parameter ranges
        num_params_r__ = 0U;
        param_ranges_i__.clear();
        num_params_r__ += t_max*k;
    }

    void set_copula_type(std::vector<int> copula_type_)
    {
        copula_type = copula_type_;

        num_params_r__ = 0U;
        param_ranges_i__.clear();
        num_params_r__ += t_max*k;
    }

    void set_latent_copula_type(std::vector<int> latent_copula_type_)
    {
        latent_copula_type = latent_copula_type_;

        num_params_r__ = 0U;
        param_ranges_i__.clear();
        num_params_r__ += t_max*k;
    }

    int get_eff_para(void){
        int eff_para = 0;
        return eff_para;
    }
    int get_t_max(void){
        return t_max;
    }

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

        Eigen::VectorXd v0_t(MCnum);
        Eigen::VectorXd weight_t(MCnum);

        v0_t[0] = ab[0] *0.5 + 0.5; // Tranform from [-1,1]
        weight_t[0] = w[0] * 0.5;// Tranform from [-1,1] to [0,1]

        for (unsigned i = 1; i < ab.size(); ++i){
            v0_t[2*i-1] = ab[i] *0.5 + 0.5;
            v0_t[2*i] = - ab[i] *0.5 + 0.5;
            weight_t[2*i-1] = w[i] * 0.5;
            weight_t[2*i] = w[i] * 0.5;
        }

        Eigen::VectorXd vg_t = v0_t;

        vector<double> logc_t(t_max,0.0);

        Eigen::VectorXd theta_12 = mean_vi.tail(eff_num_para);
        int count = 0;

        vector<double> theta(n_max,0.0);
        vector<double> theta2(n_max,0.0);

        vector<double> theta_latent(n_max,0.0);
        vector<double> theta2_latent(n_max,0.0);


        for (int i = 0; i < n_max; i++) {
            if (is_two_params(copula_type[i])) {
                theta[i] = theta_12(count); count++;
                theta2[i] = theta_12(count); count++;
            } else {
                if (copula_type[i] != 0) {
                    theta[i] = theta_12(count); count++;
                }
            }
        }

        for (int i = 0; i < n_max; i++)
        {
            if (is_two_params(latent_copula_type[i])) {
                theta_latent[i] = theta_12(count); count++;
                theta2_latent[i] = theta_12(count); count++;
            } else {
                if (latent_copula_type[i] != 0) {
                    theta_latent[i] = theta_12(count); count++;
                }
            }

        }

        for (int t = 0; t < t_max; t++) {


            Eigen::VectorXd logc_jt = Eigen::VectorXd::Zero(MCnum);

            for (int j = 0; j < MCnum; j++) { // Given v0, calculate u_cond(u,v0), j index for v0

                Eigen::VectorXd ucond_t = Eigen::VectorXd::Zero(n_max);
                Eigen::MatrixXd logc_ucond_jt = Eigen::MatrixXd::Zero(MCnum, k-1);

                for (int i = 0; i < n_max; i++) { // i index for u

                    ucond_t(i) = hfunc_trans(copula_type[i],  u(t,i), v0_t(j),  theta[i], theta2[i]);

                    for (int h = 0; h < MCnum; h++) { // h index for vg

                        logc_ucond_jt(h, gid[i]) += bicop_log_double(latent_copula_type[i],
                                                                    ucond_t(i), vg_t(h), theta_latent[i], theta2_latent[i] );
                    }


                }
                Eigen::VectorXd  max_logcvg = logc_ucond_jt.colwise().maxCoeff();
                for (int g = 0; g < k-1; g++){ // g index for group

                    Eigen::VectorXd exp_logcvg_minus_max = (logc_ucond_jt.col(g).array() - max_logcvg(g)).array().exp();
                    logc_jt(j) += max_logcvg(g) + log ( (exp_logcvg_minus_max.array() * weight_t.array() ).sum())   ;
                }



                for (int i = 0; i < n_max; i++) {
                    logc_jt(j) += bicop_log_double(copula_type[i], u(t,i), v0_t(j), theta[i], theta2[i] )   ;
                }



            }

            double max_logct = logc_jt.maxCoeff();
            Eigen::VectorXd exp_logc_jt_minus_max = (logc_jt.array() - max_logct).array().exp();
            logc_t[t] = max_logct + log ( (exp_logc_jt_minus_max.array() * weight_t.array() ).sum())   ;
        }

        for (auto& log_val : logc_t)
            logc += log_val;
        return logc;
    }

    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const
    {

        T__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        T__ DUMMY_VAR_ZERO__(0);
        (void) DUMMY_VAR__;  // suppress unused var warning

        static int current_statement_begin__;

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        // model parameters
        stan::io::reader<T__> in__(params_r__,params_i__);

        double UMAX = 1-1e-10;
        double UMIN = 1e-10;
        Eigen::Matrix<T__,Eigen::Dynamic,1>  v1;
        (void) v1;  // dummy to suppress unused var warning
        if (jacobian__)
            v1 = in__.vector_lub_constrain(UMIN,UMAX,t_max,lp__);
        else
            v1 = in__.vector_lub_constrain(UMIN,UMAX,t_max);

        Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic>  v2g;
        (void) v2g;  // dummy to suppress unused var warning
        if (jacobian__)
            v2g = in__.matrix_lub_constrain(UMIN,UMAX,t_max, k-1,lp__);
        else
            v2g = in__.matrix_lub_constrain(UMIN,UMAX,t_max, k-1);


        Eigen::Matrix<double,Eigen::Dynamic,1> u_col;
        Eigen::Matrix<T__,Eigen::Dynamic,1> u_cond_col;
        Eigen::Matrix<T__,Eigen::Dynamic,1> vg_col;

        vector_d theta_col = theta;
        vector_d theta2_col = theta2;
        vector_d latent_theta_col = latent_theta;
        vector_d latent_theta2_col = latent_theta2;

        // transformed parameters
        current_statement_begin__ = 23;
        Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic> u_cond(t_max,n_max);
        (void) u_cond;  // dummy to suppress unused var warning
        stan::math::initialize(u_cond, DUMMY_VAR__);
        stan::math::fill(u_cond,DUMMY_VAR__);



        // model body
        try
        {
            current_statement_begin__ = 12;
            int ibase = 0;

            for (int i = 0; i < n_max; i++)
            {
                current_statement_begin__ = 14;
                u_col = u.col(i);
                //VectorXd::Map(&u_col[0], t_max) = u.col(i);

                bicop_log_latent<propto__,jacobian__,double,T__,double,T__>(i, copula_type, u_col, v1, theta_col, theta2_col, lp__, lp_accum__, in__);

            } // end for

            // Transform the u copula data using Hfunc2 = u_cond
            // u_cond = hfunc_trans(copula_type, u, v1, theta, theta2 );

            for ( int i = 1; i <= n_max; i++)
            {
                int family = copula_type[i-1];
                for ( int t = 1; t <= t_max; t++){
                    stan::math::assign(get_base1_lhs(u_cond,t,i,"u_cond",1),
                                       hfunc_trans( family,
                                                    get_base1(u,t,i,"u",1),
                                                    get_base1(v1,t,"v1",1),
                                                    theta[i-1],
                                                    theta2[i-1] ));
                }
            }
            current_statement_begin__ = 13;
            ibase = 0;
            for (int i = 0; i < n_max; i++)
            {
                u_cond_col = u_cond.col(i);
                //VectorXd::Map(&u_col[0], t_max) = u.col(i);

                vg_col = v2g.col(gid[i]);
                // VectorXd::Map(&vg_col[0], t_max) = v2g.col(i);

                bicop_log_latent<propto__,jacobian__,T__,T__,double,T__>(i, latent_copula_type, u_cond_col, vg_col, latent_theta_col, latent_theta2_col, lp__, lp_accum__, in__);

            } // End for

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
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const
    {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(t_max*k);
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
        vector_d v = in__.vector_lub_constrain(UMIN,UMAX,t_max*k);

        for (int k_0__ = 0; k_0__ < t_max*k; ++k_0__)
        {
            vars__.push_back(v[k_0__]);
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
        return "bifcopulaLatent";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const
    {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= t_max*k; ++k_0__)
        {
            param_name_stream__.str(std::string());
            param_name_stream__ << "v" << '.' << k_0__;
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
        for (int k_0__ = 1; k_0__ <= t_max*k; ++k_0__)
        {
            param_name_stream__.str(std::string());
            param_name_stream__ << "v" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
    }

}; // model

} // namespace

//typedef vifcopula::bifcopulaLatent stan_model;
#endif // VIFCOPULA_BI_FACTOR_COP_LATENT_HPP

