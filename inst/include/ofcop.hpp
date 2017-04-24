#ifndef VIFCOPULA_OFCOP_HPP
#define VIFCOPULA_OFCOP_HPP

#include <omp.h>
#include <onefcopula.hpp>
#include <bicopula.hpp>
#include <stan/math.hpp>
#include <advi_mod.hpp>
#include <stan/interface_callbacks/writer/stream_writer.hpp>
#include <service/write_vb.hpp>
#include <stan/services/optimize/do_bfgs_optimize.hpp>
#include <stan/optimization/bfgs.hpp>


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


void ofcop (const matrix_d& u,
            std::vector<int>& copula_type_vec,
            int t_max,
            int n_max,
            rng_t& base_rng,
            int iter,
            int n_monte_carlo_grad,
            int n_monte_carlo_elbo,
            int eval_elbo,
            bool adapt_bool,
            double adapt_val,
            int adapt_iterations,
            double tol_rel_obj,
            int max_iterations,
            bool copselect,
            int core,
            vector_d& mean_iv, matrix_d& sample_iv, std::vector<int>& cop_vec_new)
{

    // Initiate model
    vector<double> v_temp(t_max);
    vector<double> u_temp(t_max);

    bicopula biuv(0,u_temp,v_temp,t_max,base_rng);


    int k_max = 0;
    std::vector<int> gid(n_max);
    std::fill(gid.begin(), gid.end(), 0);
    onefcopula layer_n1(u,gid,copula_type_vec,t_max, n_max, k_max, base_rng);


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
    //Could be change to Rcout in rstan
    std::stringstream out_message_writer;
    stan::interface_callbacks::writer::stream_writer message_writer(std::cout);

    std::stringstream out_parameter_writer;
    stan::interface_callbacks::writer::stream_writer parameter_writer(out_parameter_writer);

    std::stringstream out_diagnostic_writer;
    stan::interface_callbacks::writer::stream_writer diagnostic_writer(out_diagnostic_writer);


    advi_cop.run(adapt_val, adapt_bool, adapt_iterations, tol_rel_obj, 2e4,
                 message_writer, parameter_writer, diagnostic_writer);

    int max_param = layer_n1.num_params_r();
    sample_iv.resize(iter,max_param);
    mean_iv.resize(max_param);

    write_vb(out_parameter_writer, mean_iv, sample_iv);
    out_parameter_writer.clear(); // Clear state flags.


    if (copselect){
        bool keepfindcop = true;

        while (keepfindcop){
            std::cout << " Copula selection " << std::endl;
            //v_temp = mean_iv.head(t_max);
            VectorXd::Map(&v_temp[0], t_max) = mean_iv.head(t_max);

            std::vector<double> params_r(1);
            params_r[0] = 1;
            std::vector<int> params_i(0);
            bool save_iterations = false;
            int refresh = 0;
            int return_code;

            for (int j = 0; j < n_max; j++){
                //u_temp = u.col(j);
                VectorXd::Map(&u_temp[0], t_max) = u.col(j);

                biuv.reset(copula_type_vec[j], u_temp,v_temp);
                if (biuv.check_Ind()){
                    biuv.set_copula_type(0);
                    cop_vec_new[j]  = 0;

                } else {
                    const int cop_seq_size = 5;
                    int cop_seq[cop_seq_size] = {1, 3, 4, 5, 6};
                    double log_cop[cop_seq_size] = {0, 0, 0, 0, 0};
                    double AIC[cop_seq_size] = {0, 0, 0, 0, 0};
                    double BIC[cop_seq_size] = {0, 0, 0, 0, 0};
                    double lpmax = std::numeric_limits<double>::min();
                    int imax=0;
                    for (int i = 0; i < cop_seq_size; i++) {
                        biuv.set_copula_type(cop_seq[i]);
                        std::stringstream out;
                        Optimizer_BFGS bfgs(biuv, params_r, params_i, &out);
                        double lp = 0;
                        int ret = 0;
                        while (ret == 0) {
                            ret = bfgs.step();
                        }
                        lp = bfgs.logp();
                        log_cop[i] = lp;
                        AIC[i] = -2 * lp + 2 * 1;
                        BIC[i] = -2 * lp + log(t_max) * 1;
                        if (lp > lpmax){
                            lpmax = lp;
                            imax = i;
                            // std::vector<double> get_param;
                            // bfgs.params_r(get_param);
                            //cont_params[t_max+i] = get_param[0];
                        }
                    }
                    cop_vec_new[j] = cop_seq[imax];
                    cont_params[t_max+j] = 0;
                }


            }
            if (cop_vec_new != copula_type_vec){
                copula_type_vec = cop_vec_new;
                layer_n1.set_copula_type(copula_type_vec);
                advi_cop.run(adapt_val, adapt_bool, adapt_iterations, tol_rel_obj, 2e4,
                             message_writer, parameter_writer, diagnostic_writer);
            } else {
                keepfindcop = false;
            }

        }

    }

    max_param = layer_n1.num_params_r();
    sample_iv(iter,max_param);
    mean_iv(max_param);
    write_vb(out_parameter_writer, mean_iv, sample_iv);
    out_parameter_writer.clear(); // Clear state flags.





}
#endif // VIFCOPULA_OFCOP_HPP
