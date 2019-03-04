#ifndef VIFCOPULA_FVCOP_HPP
#define VIFCOPULA_FVCOP_HPP

#include <omp.h>
#include <fvcopula_stanc.hpp>
#include <fvcopulaLatent_stanc.hpp>
#include <stan/math.hpp>
#include <advi_mod.hpp>
#include <stan/callbacks/stream_writer.hpp>
#include <stan/services/optimize/bfgs.hpp>
#include <stan/optimization/bfgs.hpp>
#include <service/bicop_select.hpp>
#include <service/bicop_select_latent.hpp>
#include <transform/hfunc_stan.hpp>
#include <extra/Kruskal_MST.hpp>

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
typedef vifcopula::fvcopula fvcopula;


class fvcop
{
private:
    matrix_d u;
    vector<int> gid;             // Group of copula
    std::vector<int> copula_type;
    std::vector<int> vine_copula_type;
    matrix_int edges;

    int t_max;
    int n_max;
    int k;
    rng_t base_rng;


public:
    ~fvcop() { }

    void set_u(const matrix_d& u_)
        {
            u = u_;
        }
    void set_copula_type(std::vector<int>& copula_type_)
        {
            copula_type = copula_type_;
        }
    void set_vine_copula_type(std::vector<int>& vine_copula_type_,
                              matrix_int& edges_)
        {
            vine_copula_type = vine_copula_type_;
            edges = edges_;
        }

    fvcop (const matrix_d& u_,
            const std::vector<int>& gid_,
            std::vector<int>& copula_type_,
            std::vector<int>& vine_copula_type_,
            matrix_int& edges_,
            int t_max_,
            int n_max_,
            int k_,
            rng_t& base_rng_) :
        u(u_), gid(gid_), copula_type(copula_type_), vine_copula_type(vine_copula_type_),
        edges(edges_),
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
                    std::vector<int>& vine_cop_new,
                    matrix_int& edges_new,
                    std::vector<double>& ELBO,
                    int& count_select)
        {

            // Initiate model
            vector<double> v_temp(t_max);
            vector<double> u1_temp(t_max);
            vector<double> u2_temp(t_max);

            double ELBO_max = std::numeric_limits<double>::min();

            bicopula biuv(0,u1_temp,u2_temp,t_max,base_rng);


            fvcopula Objfvcop(u,gid,copula_type,vine_copula_type,edges,
                               t_max, n_max, k, base_rng);


            // Dummy input
            Eigen::VectorXd cont_params = Eigen::VectorXd::Zero(Objfvcop.num_params_r());

            // ADVI
            stan::variational::advi_mod<fvcopula, stan::variational::normal_meanfield, rng_t> advi_cop(Objfvcop,
                                                                                                        cont_params,
                                                                                                        base_rng,
                                                                                                        n_monte_carlo_grad,
                                                                                                        n_monte_carlo_elbo,
                                                                                                        eval_elbo,
                                                                                                        iter);
            int max_param = Objfvcop.num_params_r();
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

                    //v_temp = mean_vi.head(t_max);
                    VectorXd::Map(&v_temp[0], t_max) = mean_vi.head(t_max);

                    for (int j = 0; j < n_max; j++){
                        //u1_temp = u.col(j);
                        VectorXd::Map(&u1_temp[0], t_max) = u.col(j);
                        std::vector<double> params_out(2);
                        cop_new[j] = bicop_select(u1_temp, v_temp, t_max, params_out, base_rng);

                        for (int t = 0; t < t_max; t++){
                                u_cond(t,j) = hfunc_trans(cop_new[j],  u1_temp[t], v_temp[t], params_out[0], params_out[1]);
                            }
                        }

                    std::cout << " cop_new " << std::endl;
                    PRINT_ELEMENTS(cop_new);

                    // // Search for spaning tree here
                    matrix_d tau = MatrixXd::Zero(n_max,n_max);
                    int count_edges = 0;

                    for (int i = 0; i < n_max-1; i++){

                        VectorXd::Map(&u1_temp[0], t_max) = u_cond.col(i);

                        for (int j = i+1; j < n_max; j++){
                            if (gid[i] == gid[j]){

                                VectorXd::Map(&u2_temp[0], t_max) = u_cond.col(j);
                                bicopula biuv(1,u1_temp,u2_temp,t_max,base_rng);
                                tau(i,j) = - std::abs(biuv.kendall());
                                if (biuv.check_Ind()){
                                    tau(i,j) = 0;
                                } else{
                                // tau(j,i) = tau(i,j);
                                count_edges ++;
                                }
                            }
                        }
                    }
                    std::cout << " count_edges " << count_edges << std::endl;
                    //std::cout << " tau " << tau << std::endl;

                    matrix_int full_edges(count_edges,2);
                    matrix_int spaning_edges(n_max,2);
                    vector_d weights(count_edges);
                    int edges_id = 0;
                    for (int i = 0; i < n_max-1; i++){
                        for (int j = i+1; j < n_max; j++){
                            if (tau(i,j) != 0) {
                                full_edges(edges_id,0) = i;
                                full_edges(edges_id,1) = j;
                                weights(edges_id) = tau(i,j);
                                edges_id ++;
                            }
                        }
                    }
                    //std::cout << " full_edges " << full_edges << std::endl;
                    //std::cout << " weights " << weights << std::endl;

                    int n_group = *std::max_element(gid.begin(), gid.end()) + 1;
                    KruskalSTree(n_max, count_edges, full_edges, weights, n_group, spaning_edges);
                    vine_cop_new.resize(spaning_edges.rows());
                    // std::cout << " spaning_edges " << spaning_edges << std::endl;

                    // Bicopselect vine

                    for (int j = 0; j < spaning_edges.rows(); j++){
                        //u1_temp = u_cond.col(j);
                        VectorXd::Map(&u1_temp[0], t_max) = u_cond.col(spaning_edges(j,0));
                        VectorXd::Map(&u2_temp[0], t_max) = u_cond.col(spaning_edges(j,1));
                        std::vector<double> params_out(2);
                        vine_cop_new[j] = bicop_select_latent(u1_temp, u2_temp, t_max, params_out, base_rng);
                        }

                    std::cout << " vine_cop_new " << std::endl;
                    PRINT_ELEMENTS(vine_cop_new);



                    // TODO add check edges

                    if ((cop_new != copula_type) || (vine_cop_new != vine_copula_type) ){
                        copula_type = cop_new;
                        vine_copula_type = vine_cop_new;
                        edges_new = spaning_edges;

                        Objfvcop.set_copula_type(copula_type);
                        Objfvcop.set_vine_copula_type(vine_copula_type, edges_new);

                        max_param = Objfvcop.num_params_r();
                        cont_params = Eigen::VectorXd::Zero(max_param);

                        advi_cop.run(adapt_val, adapt_bool, adapt_iterations, tol_rel_obj, 2e4,
                                     message_writer, parameter_writer, diagnostic_writer, vi_store);
                        count_select++;

                    } else {
                        keepfindcop = false;
                    }
                    if (count_select == max_select) keepfindcop = false;

                }   // end while
            }

            stan::variational::normal_meanfield vi_save(vi_store.mu_, vi_store.omega_);
            ELBO[0] = advi_cop.calc_ELBO(vi_save, message_writer);

            max_param = Objfvcop.num_params_r();
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

        fvcopula Objfvcop(u,gid,copula_type,vine_copula_type,edges,
                            t_max, n_max, k, base_rng);

        //stan::callbacks::interrupt interrupt;
        stan::test::unit::instrumented_interrupt interrupt;

        //stan::test::unit::instrumented_logger logger;
        // stan::callbacks::stream_logger logger(debug, info, warn, error, fatal);
        stan::callbacks::stream_logger logger(std::cout,std::cout,std::cout,std::cout,std::cout);


        stan::test::unit::instrumented_writer init_writer, sample_writer, diagnostic_writer;

        std::vector<int> disc_vector;
        std::vector<double> cont_vector(Objfvcop.num_params_r() ,0); // initial value at 0

        Eigen::VectorXd inv_metric;
        stan::io::dump dmp = stan::services::util::create_unit_e_diag_inv_metric(Objfvcop.num_params_r());
        stan::io::var_context& unit_e_metric = dmp;

        try {
            inv_metric =
                stan::services::util::read_diag_inv_metric(unit_e_metric, Objfvcop.num_params_r(),
                                                           logger);
            stan::services::util::validate_diag_inv_metric(inv_metric, logger);
        } catch (const std::domain_error& e) {
        }

        stan::mcmc::adapt_diag_e_nuts<fvcopula, boost::ecuyer1988>
            sampler(Objfvcop, base_rng);

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

        stan::services::util::run_adaptive_sampler(sampler, Objfvcop, cont_vector, num_warmup,
                                                   num_samples, num_thin, refresh, save_warmup,
                                                   base_rng, interrupt, logger,
                                                   sample_writer, diagnostic_writer);


        parameter_names = sample_writer.vector_string_values();
        parameter_values = sample_writer.vector_double_values();

    } // function

    }; // class






class fvcopLatent
{
private:
    matrix_d u;
    vector_d theta;      //parameter1
    vector_d theta2;     //parameter2
    vector_d vine_theta;      //vine_theta1
    vector_d vine_theta2;     //vine_theta2
    vector<int> gid;             // Group of copula
    std::vector<int> copula_type;
    std::vector<int> vine_copula_type;
    matrix_int edges;
    int t_max;
    int n_max;
    int k;
    rng_t base_rng;


public:
    ~fvcopLatent() { }

    void set_u(const matrix_d& u_)
    {
        u = u_;
    }
    void set_copula_type(std::vector<int>& copula_type_)
    {
        copula_type = copula_type_;
    }
     void set_vine_copula_type(std::vector<int>& vine_copula_type_,
                          matrix_int& edges_)
     {
         vine_copula_type = vine_copula_type_;
         edges = edges_;
     }

    fvcopLatent (const matrix_d& u_,
                  const vector_d& theta_,
                  const vector_d& theta2_,
                  const vector_d& vine_theta_,
                  const vector_d& vine_theta2_,
            const std::vector<int>& gid_,
            std::vector<int>& copula_type_,
            std::vector<int>& vine_copula_type_,
            matrix_int& edges_,
            int t_max_,
            int n_max_,
            int k_,
            rng_t& base_rng_) :
    u(u_), theta(theta_), theta2(theta2_),
    vine_theta(vine_theta_), vine_theta2(vine_theta2_),
    gid(gid_), copula_type(copula_type_), vine_copula_type(vine_copula_type_),edges(edges_),
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
                std::vector<int>& vine_cop_new,
                matrix_int& edges_new,
                std::vector<double>& ELBO,
                int& count_select)
    {

        // Initiate model
        vector<double> v_temp(t_max);
        vector<double> u1_temp(t_max);
        vector<double> u2_temp(t_max);

        double ELBO_max = std::numeric_limits<double>::min();

        bicopula biuv(0,u1_temp,u2_temp,t_max,base_rng);


        fvcopulaLatent ObjfvcopLatent(u,theta, theta2,vine_theta, vine_theta2,
                                        gid,copula_type,vine_copula_type, edges,
                            t_max, n_max, k, base_rng);


        // Dummy input
        Eigen::VectorXd cont_params = Eigen::VectorXd::Zero(ObjfvcopLatent.num_params_r());

        // ADVI
        stan::variational::advi_mod<fvcopulaLatent, stan::variational::normal_meanfield, rng_t> advi_cop(ObjfvcopLatent,
                                                                                                    cont_params,
                                                                                                    base_rng,
                                                                                                    n_monte_carlo_grad,
                                                                                                    n_monte_carlo_elbo,
                                                                                                    eval_elbo,
                                                                                                    iter);
        int max_param = ObjfvcopLatent.num_params_r();
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

        max_param = ObjfvcopLatent.num_params_r();
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

        fvcopulaLatent ObjfvcopLatent(u,theta, theta2,vine_theta, vine_theta2,
                            gid,copula_type,vine_copula_type, edges,
                            t_max, n_max, k, base_rng);

        //stan::callbacks::interrupt interrupt;
        stan::test::unit::instrumented_interrupt interrupt;

        //stan::test::unit::instrumented_logger logger;
        // stan::callbacks::stream_logger logger(debug, info, warn, error, fatal);
        stan::callbacks::stream_logger logger(std::cout,std::cout,std::cout,std::cout,std::cout);


        stan::test::unit::instrumented_writer init_writer, sample_writer, diagnostic_writer;

        std::vector<int> disc_vector;
        std::vector<double> cont_vector(ObjfvcopLatent.num_params_r() ,0); // initial value at 0

        Eigen::VectorXd inv_metric;
        stan::io::dump dmp = stan::services::util::create_unit_e_diag_inv_metric(ObjfvcopLatent.num_params_r());
        stan::io::var_context& unit_e_metric = dmp;

        try {
            inv_metric =
                stan::services::util::read_diag_inv_metric(unit_e_metric, ObjfvcopLatent.num_params_r(),
                                                           logger);
            stan::services::util::validate_diag_inv_metric(inv_metric, logger);
        } catch (const std::domain_error& e) {
        }

        stan::mcmc::adapt_diag_e_nuts<fvcopulaLatent, boost::ecuyer1988>
            sampler(ObjfvcopLatent, base_rng);

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

        stan::services::util::run_adaptive_sampler(sampler, ObjfvcopLatent, cont_vector, num_warmup,
                                                   num_samples, num_thin, refresh, save_warmup,
                                                   base_rng, interrupt, logger,
                                                   sample_writer, diagnostic_writer);


        parameter_names = sample_writer.vector_string_values();
        parameter_values = sample_writer.vector_double_values();

    } // function

    }; // class

    }// namespace
#endif // VIFCOPULA_FVCOP_HPP
