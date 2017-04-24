#ifndef VIFCOPULA_VIBIFCOP_CPP
#define VIFCOPULA_VIBIFCOP_CPP

#include <Rcpp.h>
#include <RcppEigen.h>
#include <omp.h>
#include <ctime>
#include <onefcopula.hpp>
#include <bicopula.hpp>
#include <stan/math.hpp>
#include <advi_mod.hpp>
#include <stan/interface_callbacks/writer/stream_writer.hpp>
#include <service/write_vb.hpp>
#include <stan/services/optimize/do_bfgs_optimize.hpp>
#include <stan/optimization/bfgs.hpp>


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(StanHeaders)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(openmp)]]


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

typedef boost::ecuyer1988 rng_t;
typedef vifcopula::onefcopula onefcopula;
typedef vifcopula::bicopula bicopula;
typedef stan::optimization::BFGSLineSearch<bicopula,stan::optimization::BFGSUpdate_HInv<> > Optimizer_BFGS;

//' Variational inference for factor copula models
//'
//' \code{vifcop} returns variational estimations.
//'
//'
//' @param
//' @param
//' @return
//' @examples
//' vifcop(data, init, other)
//'
//' \dontrun{
//' vifcop(data, init, other)
//' }
//' @export
// [[Rcpp::export]]
List vibifcop(SEXP data_, SEXP init_, SEXP other_){
    BEGIN_RCPP
    static const char* function("vifcop");

    // Data input
    Rcpp::List data(data_);
        int t_max  = as<int>(data["t_max"]);
        int n_max  = as<int>(data["n_max"]);
        int k_max  = as<int>(data["k_max"]);
        matrix_d u = Rcpp::as<matrix_d>(data["u"]);
        std::vector<int> gid = Rcpp::as<std::vector<int> >(data["gid"]);
        int structfactor = as<int>(data["structfactor"]);

        stan::math::check_positive_finite(function, "Period", t_max);
        stan::math::check_positive_finite(function, "Number of variables", n_max);
        stan::math::check_positive_finite(function, "Number of latents", k_max);
        stan::math::check_positive_finite(function, "factor = 1; bifactor = 2; nestfactor = 3;", structfactor);
        //stan::math::equal(function, "Number of matrix rows",u.rows(), t_max);
        //stan::math::equal(function, "Number of matrix cols",u.cols(), n_max);
        stan::math::check_consistent_size(function, "Number of matrix columns",gid, n_max);
        stan::math::check_bounded(function, "Matrix ranges in unit space", u,0,1);
        Rcpp::Rcout << " Data input :" << " Checked" << std::endl;

    // Set other variables
    Rcpp::List other(other_);
        // Set seed
        int seed  = as<int>(other["seed"]);
        rng_t base_rng(seed);

        int core  = as<int>(other["core"]);
        // Set parallel
        #ifdef _OPENMP
            omp_set_num_threads(core);
        #endif


        int iter  = as<int>(other["iter"]); // Number of iterations after converge
        int n_monte_carlo_grad  = as<int>(other["n_monte_carlo_grad"]); //number of samples for gradient computation
        int n_monte_carlo_elbo  = as<int>(other["n_monte_carlo_elbo"]); //number of samples for ELBO computation
        int eval_elbo  = as<int>(other["eval_elbo"]);      //evaluate ELBO at every "eval_elbo" iters
        bool adapt_bool  = as<bool>(other["adapt_bool"]);      // Using adaptation
        double adapt_val  = as<double>(other["adapt_val"]);      // adaptation value
        int adapt_iterations  = as<int>(other["adapt_iterations"]);      // number of iterations for eta adaptation
        double tol_rel_obj = as<double>(other["tol_rel_obj"]);      // relative tolerance parameter for convergence
        int max_iterations = 2e4;      // max number of iterations to run algorithm
        bool copselect  = as<bool>(other["copselect"]);      // Automated copula selection

        Rcpp::Rcout << " Core : " << core << std::endl;
        Rcpp::Rcout << " General setting :" << " Checked" << std::endl;



    // Init hyperparams
    Rcpp::List init(init_);
        matrix_d v = Rcpp::as<matrix_d>(init["v"]);
        matrix_int copula_type = Rcpp::as<matrix_int>(init["copula_type"]);
        matrix_d par = Rcpp::as<matrix_d>(init["par"]);
        Rcpp::Rcout << " Init hyperparams :" << " Checked" << std::endl;


    // Timing variables
    clock_t start = clock();
    clock_t end;

    // Initiate model
    vector<double> v_temp(t_max);
    vector<double> u_temp(t_max);

    bicopula biuv(0,u_temp,v_temp,t_max,base_rng);
    vector_d par_temp = par.col(0);

    std::vector<int> copula_type_vec(n_max);
    std::vector<int> cop_vec_new(n_max);
    for (int i = 0; i < n_max; i++)
        copula_type_vec[i]= copula_type(i,0);
    k_max = 0;
    onefcopula layer_n1(u,gid,copula_type_vec,t_max, n_max, k_max, base_rng);


    // Dummy input
    Eigen::VectorXd cont_params = Eigen::VectorXd::Zero(layer_n1.num_params_r());
    // Eigen::VectorXd cont_params(t_max + par_temp.rows());
    // vector_d v_init = v.col(0);
    // cont_params << v_init, par_temp;


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
    matrix_d sample_iv(iter,max_param);
    vector_d mean_iv(max_param);
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




    Rcpp::List holder = List::create(Rcpp::Named("mean_iv") = mean_iv,
                                     Rcpp::Named("sample_iv") = sample_iv,
                                     Rcpp::Named("cop_vec_new") = cop_vec_new
    );

    end = clock();
    double delta_t = static_cast<double>(end - start) / CLOCKS_PER_SEC;

    std::cout << "It took " << delta_t << " seconds.\n"  <<  std::endl;

    // List MCMCout            = List::create(Rcpp::Named("posterior") = out_parameter_writer
    //                                            );

    return holder;
    PutRNGstate();

    END_RCPP
}
#endif // VIFCOPULA_VIBIFCOP_CPP
