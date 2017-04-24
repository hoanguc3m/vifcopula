#ifndef VIFCOPULA_VIOFCOP_CPP
#define VIFCOPULA_VIOFCOP_CPP

#include <Rcpp.h>
#include <RcppEigen.h>
#include <omp.h>
#include <ctime>
#include <stan/math.hpp>
#include <ofcop.hpp>

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
//' \code{viofcop} returns variational estimations.
//'
//'
//' @param
//' @param
//' @return
//' @examples
//' viofcop(data, init, other)
//'
//' \dontrun{
//' viofcop(data, init, other)
//' }
//' @export
// [[Rcpp::export]]
List viofcop(SEXP data_, SEXP init_, SEXP other_){
    BEGIN_RCPP
    static const char* function("viofcop");

    // Data input
    Rcpp::List data(data_);
        int t_max  = as<int>(data["t_max"]);
        int n_max  = as<int>(data["n_max"]);
        //  int k_max  = as<int>(data["k_max"]);
        matrix_d u = Rcpp::as<matrix_d>(data["u"]);
        //  std::vector<int> gid = Rcpp::as<std::vector<int> >(data["gid"]);
        //  int structfactor = as<int>(data["structfactor"]);

        stan::math::check_positive_finite(function, "Period", t_max);
        stan::math::check_positive_finite(function, "Number of variables", n_max);
        //  stan::math::check_positive_finite(function, "Number of latents", k_max);
        //  stan::math::check_positive_finite(function, "factor = 1; bifactor = 2; nestfactor = 3;", structfactor);
        //  stan::math::equal(function, "Number of matrix rows",u.rows(), t_max);
        //  stan::math::equal(function, "Number of matrix cols",u.cols(), n_max);
        //  stan::math::check_consistent_size(function, "Number of matrix columns",gid, n_max);
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

    std::vector<int> copula_type_vec(n_max);
    std::vector<int> cop_vec_new(n_max);
    for (int i = 0; i < n_max; i++)
        copula_type_vec[i]= copula_type(i,0);
    matrix_d sample_iv(iter,n_max);
    vector_d mean_iv(n_max);

    ofcop(u, copula_type_vec, t_max, n_max, base_rng,
          iter, n_monte_carlo_grad, n_monte_carlo_elbo, eval_elbo,
          adapt_bool, adapt_val, adapt_iterations, tol_rel_obj, max_iterations,
          copselect, core,
          mean_iv, sample_iv, cop_vec_new);



    Rcpp::List holder = List::create(Rcpp::Named("mean_iv") = mean_iv,
                                     Rcpp::Named("sample_iv") = sample_iv,
                                     Rcpp::Named("cop_vec_new") = cop_vec_new
    );

    end = clock();
    double delta_t = static_cast<double>(end - start) / CLOCKS_PER_SEC;

    std::cout << "It took " << delta_t << " seconds.\n"  <<  std::endl;

    return holder;
    PutRNGstate();

    END_RCPP
}
#endif // VIFCOPULA_VIOFCOP_CPP
