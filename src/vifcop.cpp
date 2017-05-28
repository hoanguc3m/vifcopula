#ifndef VIFCOPULA_VIFCOP_CPP
#define VIFCOPULA_VIFCOP_CPP

#include <Rcpp.h>
#include <RcppEigen.h>
#include <omp.h>
#include <ctime>
#include <stan/math.hpp>
#include <ofcop.hpp>
#include <nestfcop.hpp>
#include <transform/hfunc.hpp>

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

void save_vi( std::vector<string>& model_pars,
              std::vector<double>& mean_iv_save,
              matrix_d& sample_iv_save,
              const vector_d& mean_iv,
              const matrix_d& sample_iv,
              matrix_int& copula_type,
              const std::vector<int>& copula_type_vec,
              std::vector<int>& cop_vec_new,
              matrix_int& latent_copula_type,
              const std::vector<int>& latent_copula_type_vec,
              std::vector<int>& latent_cop_vec_new,
              int t_max, int n_max,
              int k, int iter, int structfactor,
              bool copselect)
{
    // Save for each tree layer
    std::string v_name;
    std::string theta_name;

    if (structfactor == 1)
    {
        v_name = "v";
        theta_name = "theta";
    }
    else
    {
        v_name = "vg";
        v_name = "thetag";
    }

    for (int i = 0; i < t_max; i++)
    {
        model_pars.push_back("v" + std::to_string(k+1) + "." + std::to_string(i+1));
        mean_iv_save.push_back(mean_iv[i]);
    }
    for (int i = t_max; i < length(mean_iv); i++)
    {
        model_pars.push_back("theta" + std::to_string(k+1) + "." + std::to_string(i+1-t_max));
        mean_iv_save.push_back(mean_iv[i]);
    }
    sample_iv_save.conservativeResize(NoChange, sample_iv_save.cols()+sample_iv.cols());
    sample_iv_save.block(0,sample_iv_save.cols()-sample_iv.cols(),iter, sample_iv.cols()) = sample_iv;

    if (copselect)
    {
        copula_type.col(0) = VectorXi::Map(&cop_vec_new[0], n_max);
        latent_copula_type.col(0) = VectorXi::Map(&latent_cop_vec_new[0], k-1);
    }
    else
    {
        cop_vec_new = copula_type_vec;
        latent_cop_vec_new = latent_copula_type_vec;
    }
}
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
List vifcop(SEXP data_, SEXP init_, SEXP other_)
{
    BEGIN_RCPP
    static const char* function("vifcop");

    // Data input
    Rcpp::List data(data_);
    int t_max  = as<int>(data["t_max"]);
    int n_max  = as<int>(data["n_max"]);
    int k_max  = as<int>(data["k_max"]);
    matrix_d u = Rcpp::as<matrix_d>(data["u"]);
    int structfactor = as<int>(data["structfactor"]);

    std::vector<int> gid = Rcpp::as<std::vector<int> >(data["gid"]);
    // Create matrix to handle group data
    int n_group = 0;
    for( int i = 0; i < n_max; i++)
    {
        if (n_group < gid[i]) n_group = gid[i];
    }
    std::vector<std::vector<int> >  g_mat(n_group, std::vector<int>(n_max));
    std::vector<int> g_count(n_group);
    for( int i = 0; i < n_max; i++)
    {
        gid[i]--;
        g_mat[gid[i]][g_count[gid[i]]] = i;
        g_count[gid[i]]++;
    }
    for( int i = 0; i < n_group; i++)
    {
        Rcpp::Rcout << " g_count " << i << ": " <<  g_count[i] << std::endl;
    }

    stan::math::check_positive_finite(function, "Period", t_max);
    stan::math::check_positive_finite(function, "Number of variables", n_max);
    stan::math::check_positive_finite(function, "Number of latents", k_max);
    stan::math::check_positive_finite(function, "factor = 1; bifactor = 2; nestfactor = 3;", structfactor);
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
    matrix_int latent_copula_type(k_max-1,1);
    matrix_d par = Rcpp::as<matrix_d>(init["par"]);
    Rcpp::Rcout << " Init hyperparams :" << " Checked" << std::endl;
    if (structfactor == 3)
    {
        latent_copula_type = Rcpp::as<matrix_int>(init["latent_copula_type"]);
    }

    // Timing variables
    clock_t start = clock();
    clock_t end;

    std::vector<int> copula_type_vec(n_max);
    std::vector<int> cop_vec_new(n_max);
    std::vector<int> latent_copula_type_vec(k_max-1);
    std::vector<int> latent_cop_vec_new(k_max-1);

    matrix_d sample_iv(iter,n_max);
    vector_d mean_iv(n_max);

    std::vector<string> model_pars;
    std::vector<double> mean_iv_save;
    matrix_d sample_iv_save(iter,0);
    switch (structfactor)
    {
    case 1: // multifactor copula
        for (int k = 0; k < k_max; k++)
        {
            // copula_type_vec = copula_type.col(k);
            VectorXi::Map(&copula_type_vec[0], n_max) = copula_type.col(k);
            if (k > 0)
            {
                hfunc_trans(u,mean_iv,cop_vec_new);
            }
            Rcpp::Rcout << "########################################################" << std::endl;
            Rcpp::Rcout << " VI Estimating copula layer:" << k+1 << std::endl;
            Rcpp::Rcout << "########################################################" << std::endl;

            ofcop Objfcop(u, copula_type_vec, t_max, n_max, k, base_rng,
                          iter, n_monte_carlo_grad, n_monte_carlo_elbo, eval_elbo,
                          adapt_bool, adapt_val, adapt_iterations, tol_rel_obj, max_iterations,
                          copselect, core);
            Objfcop.runvi(mean_iv, sample_iv, cop_vec_new);

            save_vi(model_pars, mean_iv_save, sample_iv_save,
                    mean_iv, sample_iv,
                    copula_type, copula_type_vec, cop_vec_new,
                    latent_copula_type, latent_copula_type_vec, latent_cop_vec_new,
                    t_max, n_max, k, iter, structfactor, copselect);


        }
        break;

    // case 2: // bifactor copula
    // {
    //     int k = 0;
    //     // copula_type_vec = copula_type.col(k);
    //     VectorXi::Map(&copula_type_vec[0], n_max) = copula_type.col(0);
    //     Rcpp::Rcout << "########################################################" << std::endl;
    //     Rcpp::Rcout << " VI Estimating copula layer:" << k+1 << std::endl;
    //     Rcpp::Rcout << "########################################################" << std::endl;
    //
    //     ofcop Objfcop(u, copula_type_vec, t_max, n_max, k, base_rng,
    //                   iter, n_monte_carlo_grad, n_monte_carlo_elbo, eval_elbo,
    //                   adapt_bool, adapt_val, adapt_iterations, tol_rel_obj, max_iterations,
    //                   copselect, core);
    //     Objfcop.runvi(mean_iv, sample_iv, cop_vec_new);
    //
    //     save_vi(model_pars, mean_iv_save, sample_iv_save,
    //             mean_iv, sample_iv,
    //             copula_type, copula_type_vec, cop_vec_new,
    //             latent_copula_type, latent_copula_type_vec, latent_cop_vec_new,
    //             t_max, n_max, k, iter, structfactor, copselect);
    //
    //     hfunc_trans(u,mean_iv,cop_vec_new);
    //
    //     for (int k = 1; k <= n_group; k++)
    //     {
    //         // copula_type_vec = copula_type.col(k);
    //         VectorXi::Map(&copula_type_vec[0], n_max) = copula_type.col(k);
    //
    //         Rcpp::Rcout << "########################################################" << std::endl;
    //         Rcpp::Rcout << " VI Estimating copula 2th layer: Group " << k << " with " << g_count[k-1] << " series" << std::endl;
    //         Rcpp::Rcout << "########################################################" << std::endl;
    //
    //         matrix_d u_gid(t_max,g_count[k-1]);
    //         for (int i = 0; i < g_count[k-1]; i++)
    //         {
    //             u_gid.col(i) = u.col(g_mat[k-1][i]);
    //         }
    //         ofcop Objbicop(u_gid, copula_type_vec, t_max, g_count[k-1], k, base_rng,
    //                        iter, n_monte_carlo_grad, n_monte_carlo_elbo, eval_elbo,
    //                        adapt_bool, adapt_val, adapt_iterations, tol_rel_obj, max_iterations,
    //                        copselect, core);
    //         Objbicop.runvi(mean_iv, sample_iv, cop_vec_new);
    //
    //         save_vi(model_pars, mean_iv_save, sample_iv_save,
    //                 mean_iv, sample_iv,
    //                 copula_type, copula_type_vec, cop_vec_new,
    //                 latent_copula_type, latent_copula_type_vec, latent_cop_vec_new,
    //                 t_max, g_count[k-1], k, iter, structfactor, copselect);
    //
    //
    //     }
    // }
    // break;
    case 3: // nest factor copula
    {
        int k = k_max;
        // copula_type_vec = copula_type.col(k);
        VectorXi::Map(&copula_type_vec[0], n_max) = copula_type.col(0);
        VectorXi::Map(&latent_copula_type_vec[0], k-1) = latent_copula_type.col(0);

        Rcpp::Rcout << "########################################################" << std::endl;
        Rcpp::Rcout << " VI Estimating nested factor copula" << std::endl;
        Rcpp::Rcout << "########################################################" << std::endl;

        nestfcop Objnestfcop(u, gid, copula_type_vec,latent_copula_type_vec, t_max, n_max, k, base_rng,
                             iter, n_monte_carlo_grad, n_monte_carlo_elbo, eval_elbo,
                             adapt_bool, adapt_val, adapt_iterations, tol_rel_obj, max_iterations,
                             copselect, core);
        Objnestfcop.runvi(mean_iv, sample_iv, cop_vec_new, latent_cop_vec_new);

        save_vi(model_pars, mean_iv_save, sample_iv_save,
                mean_iv, sample_iv,
                copula_type, copula_type_vec, cop_vec_new,
                latent_copula_type, latent_copula_type_vec, latent_cop_vec_new,
                t_max, n_max, k, iter, structfactor, copselect);
    }
    break;
    } // end switch



    Rcpp::List holder = List::create(Rcpp::Named("mean_iv") = mean_iv_save,
                                     Rcpp::Named("sample_iv") = sample_iv_save,
                                     Rcpp::Named("cop_type") = copula_type,
                                     Rcpp::Named("model_pars") = model_pars,
                                     Rcpp::Named("u") = u,
                                     Rcpp::Named("gid") = gid,
                                     Rcpp::Named("latent_copula_type") = latent_copula_type
                                    );

    end = clock();
    double delta_t = static_cast<double>(end - start) / CLOCKS_PER_SEC;

    std::cout << "It took " << delta_t << " seconds.\n"  <<  std::endl;

    return holder;
    PutRNGstate();

    END_RCPP
}


#endif // VIFCOPULA_VIFCOP_CPP
