#ifndef VIFCOPULA_BIFCOP_HPP
#define VIFCOPULA_BIFCOP_HPP

#include <omp.h>
#include <bifcopula_stanc.hpp>
#include <bifcopulaLatent_stanc.hpp>
#include <stan/math.hpp>
#include <advi_mod.hpp>
#include <stan/callbacks/stream_writer.hpp>
#include <stan/services/optimize/bfgs.hpp>
#include <stan/optimization/bfgs.hpp>
#include <service/bicop_select.hpp>
#include <service/bicop_select_latent.hpp>
#include <transform/hfunc_stan.hpp>

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
typedef vifcopula::bifcopula bifcopula;


class bifcop
{
private:
    matrix_d u;
    vector<int> gid;             // Group of copula
    std::vector<int> copula_type;
    std::vector<int> latent_copula_type;
    int t_max;
    int n_max;
    int k;
    rng_t base_rng;


public:
    ~bifcop() { }

    void set_u(const matrix_d& u_)
        {
            u = u_;
        }
    void set_copula_type(std::vector<int>& copula_type_)
        {
            copula_type = copula_type_;
        }
    void set_latent_copula_type(std::vector<int>& latent_copula_type_)
        {
            latent_copula_type = latent_copula_type_;
        }

    bifcop (const matrix_d& u_,
            const std::vector<int>& gid_,
            std::vector<int>& copula_type_,
            std::vector<int>& latent_copula_type_,
            int t_max_,
            int n_max_,
            int k_,
            rng_t& base_rng_) :
        u(u_), gid(gid_), copula_type(copula_type_), latent_copula_type(latent_copula_type_),
        t_max(t_max_), n_max(n_max_), k(k_),
        base_rng(base_rng_)
        {

        }

        void runvi( int iter,
                    int n_monte_carlo_grad,
                    int n_monte_carlo_elbo,
                    int eval_elbo,
                    bool adapt_bool,
                    double adapt_val,
                    int adapt_iterations,
                    double tol_rel_obj,
                    int max_iterations,
                    bool copselect,
                    bool modelselect,
                    int max_select,
                    int core,
                    vector_d& mean_vi,
                    matrix_d& sample_vi,
                    std::vector<int>& cop_new,
                    std::vector<int>& latent_cop_new,
                    std::vector<double>& ELBO,
                    int& count_select)
        {

            // Initiate model
            vector<double> v1_temp(t_max);
            vector<double> v2g_temp(t_max);
            vector<double> u_temp(t_max);
            double ELBO_max = std::numeric_limits<double>::min();

            bicopula biuv(0,u_temp,v2g_temp,t_max,base_rng);


            bifcopula Objbifcop(u,gid,copula_type,latent_copula_type,
                               t_max, n_max, k, base_rng);


            // Dummy input
            Eigen::VectorXd cont_params = Eigen::VectorXd::Zero(Objbifcop.num_params_r());

            // ADVI
            stan::variational::advi_mod<bifcopula, stan::variational::normal_meanfield, rng_t> advi_cop(Objbifcop,
                                                                                                        cont_params,
                                                                                                        base_rng,
                                                                                                        n_monte_carlo_grad,
                                                                                                        n_monte_carlo_elbo,
                                                                                                        eval_elbo,
                                                                                                        iter);
            int max_param = Objbifcop.num_params_r();
            sample_vi.resize(iter,max_param);
            mean_vi.resize(max_param);

            //Could be change to Rcout in rstan
            std::stringstream out_message_writer;
            // stan::callbacks::stream_logger logger(debug, info, warn, error, fatal);
            stan::callbacks::stream_logger message_writer(std::cout,std::cout,std::cout,std::cout,std::cout);

            std::stringstream out_parameter_writer;
            stan::callbacks::stream_writer parameter_writer(out_parameter_writer);

            std::stringstream out_diagnostic_writer;
            stan::callbacks::stream_writer diagnostic_writer(out_diagnostic_writer);

            //stan::variational::normal_meanfield vi_save(max_param);
            vifcopula::normal_meanfield_store vi_store;

            advi_cop.run(adapt_val, adapt_bool, adapt_iterations, tol_rel_obj, 2e4,
                         message_writer, parameter_writer, diagnostic_writer, vi_store);
            out_parameter_writer.clear(); // Clear state flags.


            if (copselect){
                bool keepfindcop = true;

                while (keepfindcop){
                    std::cout << "########################################################" << std::endl;
                    std::cout << " Copula selection " << std::endl;
                    std::cout << "########################################################" << std::endl;

                    stan::variational::normal_meanfield vi_tmp(vi_store.mu_, vi_store.omega_);
                    advi_cop.get_mean(vi_tmp, mean_vi);
                    matrix_d u_cond(t_max,n_max);
                    //v1_temp = mean_vi.head(t_max);
                    matrix_d v2g = mean_vi.head(t_max*k);
                    VectorXd::Map(&v1_temp[0], t_max) = mean_vi.head(t_max);

                    v2g.resize(t_max,k);

                    for (int j = 0; j < n_max; j++){
                        //u_temp = u.col(j);
                        VectorXd::Map(&u_temp[0], t_max) = u.col(j);
                        std::vector<double> params_out(2);
                        cop_new[j] = bicop_select(u_temp, v1_temp, t_max, params_out, base_rng);

                        for (int t = 0; t < t_max; t++){
                                int cop_temp = cop_new[j];
                                if (cop_temp == 21 || cop_temp == 22 || cop_temp == 25) {
                                    cop_temp = cop_temp - 20;
                                }
                                u_cond(t,j) = hfunc_trans(cop_temp,  u_temp[t], v1_temp[t], params_out[0], params_out[1]);
                                // if ( u_cond(t,j) < 1.1e-10 ||  u_cond(t,j) > 1 - 1.1e-10 ){
                                //     Rcpp::Rcout << " u " << u_cond(t,j) << " " << cop_new[j] << " " <<
                                //         u_temp[t] << " " << v1_temp[t] << " " << params_out[0] << " " << params_out[1] << std::endl;
                                // }
                            }
                        }


                    for (int j = 0; j < n_max; j++){
                        //u_temp = u_cond.col(j);
                        VectorXd::Map(&u_temp[0], t_max) = u_cond.col(j);
                        VectorXd::Map(&v2g_temp[0], t_max) = v2g.col(gid[j]+1);
                        std::vector<double> params_out(2);
                        latent_cop_new[j] = bicop_select_latent(u_temp, v2g_temp, t_max, params_out, base_rng);
                        }


                    std::cout << " cop_new " << std::endl;
                    PRINT_ELEMENTS(cop_new);

                    std::cout << " latent_cop_new " << std::endl;
                    PRINT_ELEMENTS(latent_cop_new);



                    if ((cop_new != copula_type) || (latent_cop_new != latent_copula_type)){
                        copula_type = cop_new;
                        latent_copula_type = latent_cop_new;

                        Objbifcop.set_copula_type(copula_type);
                        Objbifcop.set_latent_copula_type(latent_copula_type);

                        max_param = Objbifcop.num_params_r();
                        cont_params = Eigen::VectorXd::Zero(max_param);

                        advi_cop.run(adapt_val, adapt_bool, adapt_iterations, tol_rel_obj, 2e4,
                                     message_writer, parameter_writer, diagnostic_writer, vi_store);
                        count_select++;

                        // stan::variational::normal_meanfield vi_save(vi_store.mu_, vi_store.omega_);
                        // ELBO[0] = advi_cop.calc_ELBO(vi_save, message_writer);
                        // if (ELBO[0] < ELBO_max){
                        //     if( abs(ELBO[0] / ELBO_max - 1) < 0.001 ) { // stop until convergence
                        //         keepfindcop = false;
                        //     } else {
                        //         ELBO_max = ELBO[0];
                        //     }
                        //
                        // } else{
                        //     ELBO_max = ELBO[0];
                        // }
                    } else {
                        keepfindcop = false;
                    }
                    if (count_select == max_select) keepfindcop = false;

                }   // end while
            }



            stan::variational::normal_meanfield vi_save(vi_store.mu_, vi_store.omega_);
            ELBO[0] = advi_cop.calc_ELBO(vi_save, message_writer);

            max_param = Objbifcop.num_params_r();
            mean_vi.resize(max_param);
            sample_vi.resize(iter,max_param);
            advi_cop.write(vi_save, mean_vi, sample_vi, ELBO, modelselect, message_writer);
            out_parameter_writer.clear(); // Clear state flags.
            std::cout << " Done ! " << std::endl;





        } // function

    void runhmc( int num_warmup,
                 int num_samples,
                 int num_thin,
                 bool save_warmup,
                 int refresh,
                 unsigned int chain,
                 double init_radius,
                 double stepsize,
                 double stepsize_jitter,
                 int max_depth,
                 double delta,
                 double gamma,
                 double kappa,
                 double t0,
                 unsigned int init_buffer,
                 unsigned int term_buffer,
                 unsigned int window,
                 std::vector<std::vector<std::string> >& parameter_names,
                 std::vector<std::vector<double> >& parameter_values
    )
    {

        bifcopula Objbifcop(u,gid,copula_type,latent_copula_type,
                            t_max, n_max, k, base_rng);

        //stan::callbacks::interrupt interrupt;
        stan::test::unit::instrumented_interrupt interrupt;

        //stan::test::unit::instrumented_logger logger;
        // stan::callbacks::stream_logger logger(debug, info, warn, error, fatal);
        stan::callbacks::stream_logger logger(std::cout,std::cout,std::cout,std::cout,std::cout);


        stan::test::unit::instrumented_writer init_writer, sample_writer, diagnostic_writer;

        std::vector<int> disc_vector;
        std::vector<double> cont_vector(Objbifcop.num_params_r() ,0); // initial value at 0

        Eigen::VectorXd inv_metric;
        stan::io::dump dmp = stan::services::util::create_unit_e_diag_inv_metric(Objbifcop.num_params_r());
        stan::io::var_context& unit_e_metric = dmp;

        try {
            inv_metric =
                stan::services::util::read_diag_inv_metric(unit_e_metric, Objbifcop.num_params_r(),
                                                           logger);
            stan::services::util::validate_diag_inv_metric(inv_metric, logger);
        } catch (const std::domain_error& e) {
        }

        stan::mcmc::adapt_diag_e_nuts<bifcopula, boost::ecuyer1988>
            sampler(Objbifcop, base_rng);

        sampler.set_metric(inv_metric);
        sampler.set_nominal_stepsize(stepsize);
        sampler.set_stepsize_jitter(stepsize_jitter);
        sampler.set_max_depth(max_depth);

        sampler.get_stepsize_adaptation().set_mu(log(10 * stepsize));
        sampler.get_stepsize_adaptation().set_delta(delta);
        sampler.get_stepsize_adaptation().set_gamma(gamma);
        sampler.get_stepsize_adaptation().set_kappa(kappa);
        sampler.get_stepsize_adaptation().set_t0(t0);

        sampler.set_window_params(num_warmup, init_buffer, term_buffer,
                                  window, logger);

        stan::services::util::run_adaptive_sampler(sampler, Objbifcop, cont_vector, num_warmup,
                                                   num_samples, num_thin, refresh, save_warmup,
                                                   base_rng, interrupt, logger,
                                                   sample_writer, diagnostic_writer);


        parameter_names = sample_writer.vector_string_values();
        parameter_values = sample_writer.vector_double_values();

    } // function

    }; // class






class bifcopLatent
{
private:
    matrix_d u;
    vector_d theta;      //parameter1
    vector_d theta2;     //parameter2
    vector_d latent_theta;      //latent_theta1
    vector_d latent_theta2;     //latent_theta2
    vector<int> gid;             // Group of copula
    std::vector<int> copula_type;
    std::vector<int> latent_copula_type;
    int t_max;
    int n_max;
    int k;
    rng_t base_rng;


public:
    ~bifcopLatent() { }

    void set_u(const matrix_d& u_)
    {
        u = u_;
    }
    void set_copula_type(std::vector<int>& copula_type_)
    {
        copula_type = copula_type_;
    }
    void set_latent_copula_type(std::vector<int>& latent_copula_type_)
    {
        latent_copula_type = latent_copula_type_;
    }

    bifcopLatent (const matrix_d& u_,
                  const vector_d& theta_,
                  const vector_d& theta2_,
                  const vector_d& latent_theta_,
                  const vector_d& latent_theta2_,
            const std::vector<int>& gid_,
            std::vector<int>& copula_type_,
            std::vector<int>& latent_copula_type_,
            int t_max_,
            int n_max_,
            int k_,
            rng_t& base_rng_) :
    u(u_), theta(theta_), theta2(theta2_),
    latent_theta(latent_theta_), latent_theta2(latent_theta2_),
    gid(gid_), copula_type(copula_type_), latent_copula_type(latent_copula_type_),
    t_max(t_max_), n_max(n_max_), k(k_),
    base_rng(base_rng_)
    {

    }

    void runvi( int iter,
                int n_monte_carlo_grad,
                int n_monte_carlo_elbo,
                int eval_elbo,
                bool adapt_bool,
                double adapt_val,
                int adapt_iterations,
                double tol_rel_obj,
                int max_iterations,
                bool copselect,
                bool modelselect,
                int max_select,
                int core,
                vector_d& mean_vi,
                matrix_d& sample_vi,
                std::vector<int>& cop_new,
                std::vector<int>& latent_cop_new,
                std::vector<double>& ELBO,
                int& count_select)
    {

        // Initiate model
        vector<double> v1_temp(t_max);
        vector<double> v2g_temp(t_max);
        vector<double> u_temp(t_max);
        double ELBO_max = std::numeric_limits<double>::min();

        bicopula biuv(0,u_temp,v2g_temp,t_max,base_rng);


        bifcopulaLatent ObjbifcopLatent(u,theta, theta2,latent_theta, latent_theta2,
                                        gid,copula_type,latent_copula_type,
                            t_max, n_max, k, base_rng);


        // Dummy input
        Eigen::VectorXd cont_params = Eigen::VectorXd::Zero(ObjbifcopLatent.num_params_r());

        // ADVI
        stan::variational::advi_mod<bifcopulaLatent, stan::variational::normal_meanfield, rng_t> advi_cop(ObjbifcopLatent,
                                                                                                    cont_params,
                                                                                                    base_rng,
                                                                                                    n_monte_carlo_grad,
                                                                                                    n_monte_carlo_elbo,
                                                                                                    eval_elbo,
                                                                                                    iter);
        int max_param = ObjbifcopLatent.num_params_r();
        sample_vi.resize(iter,max_param);
        mean_vi.resize(max_param);

        //Could be change to Rcout in rstan
        std::stringstream out_message_writer;
        // stan::callbacks::stream_logger logger(debug, info, warn, error, fatal);
        stan::callbacks::stream_logger message_writer(std::cout,std::cout,std::cout,std::cout,std::cout);

        std::stringstream out_parameter_writer;
        stan::callbacks::stream_writer parameter_writer(out_parameter_writer);

        std::stringstream out_diagnostic_writer;
        stan::callbacks::stream_writer diagnostic_writer(out_diagnostic_writer);

        //stan::variational::normal_meanfield vi_save(max_param);
        vifcopula::normal_meanfield_store vi_store;

        advi_cop.run(adapt_val, adapt_bool, adapt_iterations, tol_rel_obj, 2e4,
                     message_writer, parameter_writer, diagnostic_writer, vi_store);
        out_parameter_writer.clear(); // Clear state flags.


        stan::variational::normal_meanfield vi_save(vi_store.mu_, vi_store.omega_);
        ELBO[0] = advi_cop.calc_ELBO(vi_save, message_writer);

        max_param = ObjbifcopLatent.num_params_r();
        mean_vi.resize(max_param);
        sample_vi.resize(iter,max_param);
        modelselect = false; // no calculation for Latent infer only
        advi_cop.write(vi_save, mean_vi, sample_vi, ELBO, modelselect, message_writer);
        out_parameter_writer.clear(); // Clear state flags.
        std::cout << " Done ! " << std::endl;





    } // function

    void runhmc( int num_warmup,
                 int num_samples,
                 int num_thin,
                 bool save_warmup,
                 int refresh,
                 unsigned int chain,
                 double init_radius,
                 double stepsize,
                 double stepsize_jitter,
                 int max_depth,
                 double delta,
                 double gamma,
                 double kappa,
                 double t0,
                 unsigned int init_buffer,
                 unsigned int term_buffer,
                 unsigned int window,
                 std::vector<std::vector<std::string> >& parameter_names,
                 std::vector<std::vector<double> >& parameter_values
    )
    {

        bifcopulaLatent ObjbifcopLatent(u,theta, theta2,latent_theta, latent_theta2,
                            gid,copula_type,latent_copula_type,
                            t_max, n_max, k, base_rng);

        //stan::callbacks::interrupt interrupt;
        stan::test::unit::instrumented_interrupt interrupt;

        //stan::test::unit::instrumented_logger logger;
        // stan::callbacks::stream_logger logger(debug, info, warn, error, fatal);
        stan::callbacks::stream_logger logger(std::cout,std::cout,std::cout,std::cout,std::cout);


        stan::test::unit::instrumented_writer init_writer, sample_writer, diagnostic_writer;

        std::vector<int> disc_vector;
        std::vector<double> cont_vector(ObjbifcopLatent.num_params_r() ,0); // initial value at 0

        Eigen::VectorXd inv_metric;
        stan::io::dump dmp = stan::services::util::create_unit_e_diag_inv_metric(ObjbifcopLatent.num_params_r());
        stan::io::var_context& unit_e_metric = dmp;

        try {
            inv_metric =
                stan::services::util::read_diag_inv_metric(unit_e_metric, ObjbifcopLatent.num_params_r(),
                                                           logger);
            stan::services::util::validate_diag_inv_metric(inv_metric, logger);
        } catch (const std::domain_error& e) {
        }

        stan::mcmc::adapt_diag_e_nuts<bifcopulaLatent, boost::ecuyer1988>
            sampler(ObjbifcopLatent, base_rng);

        sampler.set_metric(inv_metric);
        sampler.set_nominal_stepsize(stepsize);
        sampler.set_stepsize_jitter(stepsize_jitter);
        sampler.set_max_depth(max_depth);

        sampler.get_stepsize_adaptation().set_mu(log(10 * stepsize));
        sampler.get_stepsize_adaptation().set_delta(delta);
        sampler.get_stepsize_adaptation().set_gamma(gamma);
        sampler.get_stepsize_adaptation().set_kappa(kappa);
        sampler.get_stepsize_adaptation().set_t0(t0);

        sampler.set_window_params(num_warmup, init_buffer, term_buffer,
                                  window, logger);

        stan::services::util::run_adaptive_sampler(sampler, ObjbifcopLatent, cont_vector, num_warmup,
                                                   num_samples, num_thin, refresh, save_warmup,
                                                   base_rng, interrupt, logger,
                                                   sample_writer, diagnostic_writer);


        parameter_names = sample_writer.vector_string_values();
        parameter_values = sample_writer.vector_double_values();

    } // function

    }; // class

    }// namespace
#endif // VIFCOPULA_BIFCOP_HPP
