#ifndef VIFCOPULA_VIFCOP_CPP
#define VIFCOPULA_VIFCOP_CPP

#include <Rcpp.h>
#include <RcppEigen.h>
#include <omp.h>
#include <ctime>
#include <one_factor_cop.cpp>
#include <rstan/rstaninc.hpp>

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
        vector_int gid = Rcpp::as<vector_int>(data["gid"]);
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
    vector_d copula_type_temp = copula_type.col(0);
    vifcopula::one_factor_cop copula_level_1(u,gid,copula_type_temp,t_max, n_max, k_max);

    // Dummy input
    //Eigen::VectorXd cont_params = Eigen::VectorXd::Zero(my_model.num_params_r());
    stan::variational::normal_meanfield cont_params(t_max + n_max);


    // ADVI
    stan::variational::advi<Model_cp, stan::variational::normal_meanfield, rng_t> test_advi(my_model,
                                                                                            cont_params,
                                                                                            base_rng,
                                                                                            10,
                                                                                            100,
                                                                                            100,
                                                                                            1);

    stan::callbacks::writer writer;

    //test_advi.run(0.01, false, 50, 1, 2e4,
    //              writer, writer, writer);


    Rcpp::Rcout << " copula LL :" << copula_level_1.log_prob() << std::endl;
    Rcpp::Rcout << " normal_meanfield LL :" << cont_params.entropy() << std::endl;





    end = clock();
    double delta_t = static_cast<double>(end - start) / CLOCKS_PER_SEC;

    std::cout << "It took " << delta_t << " seconds.\n"  <<  std::endl;

    List MCMCout            = List::create(Rcpp::Named("delta_t") = 0
                                               );
    return MCMCout;
    PutRNGstate();

    END_RCPP
}
#endif // VIFCOPULA_VIFCOP_CPP
