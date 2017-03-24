#ifndef VIFCOPULA_VIFCOP_CPP
#define VIFCOPULA_VIFCOP_CPP

#include <Rcpp.h>
#include <RcppEigen.h>
#include <omp.h>
#include <ctime>
#include <one_factor_cop.hpp>
#include <stan/math.hpp>
#include <advi_mod.hpp>
#include <stan/interface_callbacks/writer/stream_writer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


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
typedef vifcopula::one_factor_cop Model_cp;

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
List vifcop(SEXP data_, SEXP init_, SEXP other_){
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
        double adapt_val  = as<double>(other["adapt_val"]);      // Using adaptation
        int adapt_iterations  = as<int>(other["adapt_iterations"]);      // number of iterations for eta adaptation
        double tol_rel_obj = as<double>(other["tol_rel_obj"]);      // relative tolerance parameter for convergence
        int max_iterations = 2e4;      // max number of iterations to run algorithm

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
    vector_d v_temp = v.col(0);
    std::vector<int> copula_type_vec(n_max);
    for (int i = 0; i < n_max; i++)
        copula_type_vec[i]= copula_type(i,0);
    k_max = 0;
    Model_cp copula_l1(u,gid,copula_type_vec,t_max, n_max, k_max, base_rng);


    // Dummy input
    Eigen::VectorXd cont_params = Eigen::VectorXd::Zero(copula_l1.num_params_r());
    //stan::variational::normal_meanfield cont_params(t_max + n_max);


    // ADVI
    stan::variational::advi_mod<Model_cp, stan::variational::normal_meanfield, rng_t> advi_cop(copula_l1,
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

    int max_param = copula_l1.num_params_r();
    matrix_d sample_iv(iter,max_param);
    vector_d mean_iv(max_param);

    std::string token;
    //std::getline(out_parameter_writer, token);
    //std::getline(out_parameter_writer, token);

    int i,j;
    try{
        std::getline(out_parameter_writer, token, ',');
        for (j = 0; j < max_param-1;j++){
            std::getline(out_parameter_writer, token, ',');
            boost::trim(token);
            mean_iv(j) = boost::lexical_cast<double>(token);
        }
        std::getline(out_parameter_writer, token);
        boost::trim(token);
        mean_iv(max_param-1) = boost::lexical_cast<double>(token);

        for (i = 0; i < iter;i++){
            std::getline(out_parameter_writer, token, ',');
            for (j = 0; j < max_param-1;j++){
                std::getline(out_parameter_writer, token, ',');
                boost::trim(token);
                sample_iv(i,j) = boost::lexical_cast<double>(token);
            }
            std::getline(out_parameter_writer, token);
            boost::trim(token);
            sample_iv(i,max_param-1) = boost::lexical_cast<double>(token);
        }

    } catch (...){
        std::cout << "Error at i = " << i << " j = " << j << std::endl;
        std::cout << " token " << token << std::endl;
    }



    Rcpp::List holder = List::create(Rcpp::Named("mean_iv") = mean_iv,
                                     Rcpp::Named("sample_iv") = sample_iv
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
#endif // VIFCOPULA_VIFCOP_CPP
