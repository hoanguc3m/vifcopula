#ifndef VIFCOPULA_OFCOP_HPP
#define VIFCOPULA_OFCOP_HPP

#include <omp.h>
#include <onefcopula_stanc.hpp>
#include <stan/math.hpp>
#include <advi_mod.hpp>
#include <stan/callbacks/writer.hpp>
#include <stan/callbacks/stream_logger.hpp>
#include <stan/callbacks/stream_writer.hpp>
#include <stan/callbacks/interrupt.hpp>
#include <stan/services/sample/hmc_nuts_diag_e_adapt.hpp>
#include <stan/io/empty_var_context.hpp>
#include <service/instrumented_callbacks.hpp>

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


class ofcop
{
private:
    matrix_d u;
    std::vector<int> copula_type_vec;
    int t_max;
    int n_max;
    int k;
    rng_t base_rng;


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
            rng_t& base_rng_) :
    u(u_), copula_type_vec(copula_type_vec_),t_max(t_max_), n_max(n_max_), k(k_),
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
                int core,
                vector_d& mean_vi,
                matrix_d& sample_vi,
                std::vector<int>& cop_vec_new,
                std::vector<double>& ELBO,
                int& count_select)
    {

        // Initiate model
        vector<double> v_temp(t_max);
        vector<double> u_temp(t_max);
        double ELBO_max = std::numeric_limits<double>::min();


        std::vector<int> gid(n_max);
        std::fill(gid.begin(), gid.end(), 0);
        onefcopula ObjOnefcop(u,gid,copula_type_vec,t_max, n_max, k, base_rng);


        // Dummy input
        Eigen::VectorXd cont_params = Eigen::VectorXd::Zero(ObjOnefcop.num_params_r());

        // ADVI
        stan::variational::advi_mod<onefcopula, stan::variational::normal_meanfield, rng_t> advi_cop(ObjOnefcop,
                                                                                                     cont_params,
                                                                                                     base_rng,
                                                                                                     n_monte_carlo_grad,
                                                                                                     n_monte_carlo_elbo,
                                                                                                     eval_elbo,
                                                                                                     iter);
        int max_param = ObjOnefcop.num_params_r();
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

//        stan::variational::normal_meanfield vi_save(max_param);
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

                //v_temp = mean_vi.head(t_max);
                VectorXd::Map(&v_temp[0], t_max) = mean_vi.head(t_max);

                std::vector<double> params_out(2);
                matrix_d u_omp(u) ;
                int t_max_omp = t_max;
                int n_max_omp = n_max;
                rng_t base_rng_omp(0);

                // omp_set_num_threads(1);

                // #pragma omp parallel for default(none) firstprivate(u_temp,v_temp,params_out,t_max_omp, base_rng_omp) shared(n_max_omp,cop_vec_new,u_omp)
                    for (int j = 0; j < n_max_omp; j++){
                        //u_temp = u.col(j);
                        VectorXd::Map(&u_temp[0], t_max_omp) = u_omp.col(j);
                        cop_vec_new[j] = bicop_select(u_temp, v_temp, t_max_omp,params_out, base_rng_omp);
                    }

                if (cop_vec_new != copula_type_vec){
                    copula_type_vec = cop_vec_new;
                    ObjOnefcop.set_copula_type(copula_type_vec);
                    max_param = ObjOnefcop.num_params_r();
                    cont_params = Eigen::VectorXd::Zero(max_param);
                    //vi_save.resize(max_param);

                    advi_cop.run(adapt_val, adapt_bool, adapt_iterations, tol_rel_obj, 2e4,
                                 message_writer, parameter_writer, diagnostic_writer, vi_store);
                    count_select++;
                    stan::variational::normal_meanfield vi_save(vi_store.mu_, vi_store.omega_);
                    ELBO[0] = advi_cop.calc_ELBO(vi_save, message_writer);
                    if (ELBO[0] < ELBO_max){
                        if( abs(ELBO[0] / ELBO_max - 1) < 0.01 ) { // stop until convergence
                            keepfindcop = false;
                        } else {
                            ELBO_max = ELBO[0];
                        }

                    } else{
                        ELBO_max = ELBO[0];
                    }

                } else {
                    keepfindcop = false;
                }

            }

        }

        stan::variational::normal_meanfield vi_save(vi_store.mu_, vi_store.omega_);
        ELBO[0] = advi_cop.calc_ELBO(vi_save, message_writer);

        max_param = ObjOnefcop.num_params_r();
        sample_vi.resize(iter,max_param);
        mean_vi.resize(max_param);
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

        std::vector<int> gid(n_max);
        std::fill(gid.begin(), gid.end(), 0);
        onefcopula ObjOnefcop(u,gid,copula_type_vec,t_max, n_max, k, base_rng);

        //stan::callbacks::interrupt interrupt;
        stan::test::unit::instrumented_interrupt interrupt;

        //stan::test::unit::instrumented_logger logger;
        // stan::callbacks::stream_logger logger(debug, info, warn, error, fatal);
        stan::callbacks::stream_logger logger(std::cout,std::cout,std::cout,std::cout,std::cout);


        stan::test::unit::instrumented_writer init_writer, sample_writer, diagnostic_writer;

        std::vector<int> disc_vector;
        std::vector<double> cont_vector(ObjOnefcop.num_params_r() ,0); // initial value at 0

        Eigen::VectorXd inv_metric;
        stan::io::dump dmp = stan::services::util::create_unit_e_diag_inv_metric(ObjOnefcop.num_params_r());
        stan::io::var_context& unit_e_metric = dmp;

        try {
            inv_metric =
                stan::services::util::read_diag_inv_metric(unit_e_metric, ObjOnefcop.num_params_r(),
                                                           logger);
            stan::services::util::validate_diag_inv_metric(inv_metric, logger);
        } catch (const std::domain_error& e) {
        }

        stan::mcmc::adapt_diag_e_nuts<onefcopula, boost::ecuyer1988>
            sampler(ObjOnefcop, base_rng);

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

        stan::services::util::run_adaptive_sampler(sampler, ObjOnefcop, cont_vector, num_warmup,
                                                   num_samples, num_thin, refresh, save_warmup,
                                                   base_rng, interrupt, logger,
                                                   sample_writer, diagnostic_writer);


        parameter_names = sample_writer.vector_string_values();
        parameter_values = sample_writer.vector_double_values();

    } // function

}; // class

}// namespace
#endif // VIFCOPULA_OFCOP_HPP
