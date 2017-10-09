#ifndef VIFCOPULA_OFCOP_HPP
#define VIFCOPULA_OFCOP_HPP

#include <omp.h>
#include <onefcopula_stanc.hpp>
#include <stan/math.hpp>
#include <advi_mod.hpp>
#include <stan/callbacks/writer.hpp>
#include <stan/callbacks/stream_logger.hpp>
#include <service/bicop_select.hpp>
#include <service/normal_meanfield_store.hpp>

namespace vifcopula
{


using namespace Eigen;
using namespace stan::math;
using namespace vifcopula;

typedef Eigen::Matrix<int,Eigen::Dynamic,1> vector_int;
typedef Eigen::Matrix<int,1,Eigen::Dynamic> row_vector_int;
typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> matrix_int;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;

typedef boost::ecuyer1988 rng_t;
typedef vifcopula::onefcopula onefcopula;
typedef vifcopula::bicopula bicopula;
typedef stan::optimization::BFGSLineSearch<bicopula,stan::optimization::BFGSUpdate_HInv<> > Optimizer_BFGS;


class ofcop
{
private:
    matrix_d u;
    std::vector<int> copula_type_vec;
    int t_max;
    int n_max;
    int k;
    rng_t base_rng;
    int iter;
    int n_monte_carlo_grad;
    int n_monte_carlo_elbo;
    int eval_elbo;
    bool adapt_bool;
    double adapt_val;
    int adapt_iterations;
    double tol_rel_obj;
    int max_iterations;
    bool copselect;
    int core;

public:
    ~ofcop() { }

    void set_u(const matrix_d& u_)
    {
        u = u_;
    }
    void set_copula_type(std::vector<int>& copula_type_vec_)
    {
        copula_type_vec = copula_type_vec_;
    }

    ofcop (const matrix_d& u_,
            std::vector<int>& copula_type_vec_,
            int t_max_,
            int n_max_,
            int k_,
            rng_t& base_rng_,
            int iter_,
            int n_monte_carlo_grad_,
            int n_monte_carlo_elbo_,
            int eval_elbo_,
            bool adapt_bool_,
            double adapt_val_,
            int adapt_iterations_,
            double tol_rel_obj_,
            int max_iterations_,
            bool copselect_,
            int core_) :
    u(u_), copula_type_vec(copula_type_vec_),t_max(t_max_), n_max(n_max_), k(k_),
    base_rng(base_rng_), iter(iter_), n_monte_carlo_grad(n_monte_carlo_grad_),
    n_monte_carlo_elbo(n_monte_carlo_elbo_), eval_elbo(eval_elbo_),
    adapt_bool(adapt_bool_), adapt_val(adapt_val_), adapt_iterations(adapt_iterations_),
    tol_rel_obj(tol_rel_obj_), max_iterations(max_iterations_),
    copselect(copselect_), core(core_)
    {

    }
    void runvi( vector_d& mean_iv,
                matrix_d& sample_iv,
                std::vector<int>& cop_vec_new,
                double& ELBO)
    {

        // Initiate model
        vector<double> v_temp(t_max);
        vector<double> u_temp(t_max);



        std::vector<int> gid(n_max);
        std::fill(gid.begin(), gid.end(), 0);
        onefcopula layer_n1(u,gid,copula_type_vec,t_max, n_max, k, base_rng);


        // Dummy input
        Eigen::VectorXd cont_params = Eigen::VectorXd::Zero(layer_n1.num_params_r());

        // ADVI
        stan::variational::advi_mod<onefcopula, stan::variational::normal_meanfield, rng_t> advi_cop(layer_n1,
                                                                                                     cont_params,
                                                                                                     base_rng,
                                                                                                     n_monte_carlo_grad,
                                                                                                     n_monte_carlo_elbo,
                                                                                                     eval_elbo,
                                                                                                     iter);
        int max_param = layer_n1.num_params_r();
        sample_iv.resize(iter,max_param);
        mean_iv.resize(max_param);

        //Could be change to Rcout in rstan
        std::stringstream out_message_writer;
        // stan::callbacks::stream_logger logger(debug, info, warn, error, fatal);
        stan::callbacks::stream_logger message_writer(std::cout,std::cout,std::cout,std::cout,std::cout);

        std::stringstream out_parameter_writer;
        stan::callbacks::stream_writer parameter_writer(out_parameter_writer);

        std::stringstream out_diagnostic_writer;
        stan::callbacks::stream_writer diagnostic_writer(out_diagnostic_writer);

//        stan::variational::normal_meanfield vi_save(max_param);
        vifcopula::normal_meanfield_store vi_store;

        advi_cop.run(adapt_val, adapt_bool, adapt_iterations, tol_rel_obj, 2e4,
                     message_writer, parameter_writer, diagnostic_writer, vi_store);
        out_parameter_writer.clear(); // Clear state flags.


        if (copselect){
            bool keepfindcop = true;
            stan::variational::normal_meanfield vi_tmp(vi_store.mu_, vi_store.omega_);
            advi_cop.get_mean(vi_tmp, mean_iv);


            while (keepfindcop){
                std::cout << "########################################################" << std::endl;
                std::cout << " Copula selection " << std::endl;
                std::cout << "########################################################" << std::endl;
                //v_temp = mean_iv.head(t_max);
                VectorXd::Map(&v_temp[0], t_max) = mean_iv.head(t_max);


                for (int j = 0; j < n_max; j++){
                    //u_temp = u.col(j);
                    VectorXd::Map(&u_temp[0], t_max) = u.col(j);

                    std::vector<double> params_out(2);
                    cop_vec_new[j] = bicop_select(u_temp, v_temp, t_max,params_out, base_rng);
                    // cont_params[t_max+j] = 0;


                }

                // std::cout << " cop_vec_new " << std::endl;
                // PRINT_ELEMENTS(cop_vec_new);
                // keepfindcop = false;
                if (cop_vec_new != copula_type_vec){
                    copula_type_vec = cop_vec_new;
                    layer_n1.set_copula_type(copula_type_vec);
                    max_param = layer_n1.num_params_r();
                    cont_params = Eigen::VectorXd::Zero(max_param);
                    //vi_save.resize(max_param);

                    advi_cop.run(adapt_val, adapt_bool, adapt_iterations, tol_rel_obj, 2e4,
                                 message_writer, parameter_writer, diagnostic_writer, vi_store);

                } else {
                    keepfindcop = false;
                }

            }

        }

        stan::variational::normal_meanfield vi_save(vi_store.mu_, vi_store.omega_);
        ELBO = advi_cop.calc_ELBO(vi_save, message_writer);

        max_param = layer_n1.num_params_r();
        sample_iv.resize(iter,max_param);
        mean_iv.resize(max_param);
        advi_cop.write(vi_save, mean_iv, sample_iv, message_writer);


        out_parameter_writer.clear(); // Clear state flags.
        std::cout << " Done ! " << std::endl;





    } // function

}; // class

}// namespace
#endif // VIFCOPULA_OFCOP_HPP
