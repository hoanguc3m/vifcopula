////////////////////////////////////////////////////////////
// A collection of functions
////////////////////////////////////////////////////////////

#ifndef MISCELLANEOUS_HPP
#define MISCELLANEOUS_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
#include <stan/math.hpp>

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

namespace vifcopula{




//////////////////////////////////////////////////////////////////
// Change the eigen type to c++ std type ///////////////
/////////////////////////////////////////////////////////////////

void save_vi( std::vector<string>& model_pars,
              std::vector<double>& mean_vi_save,
              matrix_d& sample_vi_save,
              const vector_d& mean_vi,
              const matrix_d& sample_vi,
              std::vector<int>& copula_type,
              std::vector<int>& cop_vec_new,
              std::vector<int>& latent_copula_type,
              std::vector<int>& latent_cop_vec_new,
              int t_max, int n_max,
              int k, int iter, int structfactor,
              bool copselect)
{
    std::string v_name;
    std::string theta_name;

    // save model_pars, mean_vi_save, sample_vi_save,
    // copula_type, latent_copula_type
    if (structfactor > 10){     // Only infer v
        if (structfactor == 11 || structfactor == 14){
            v_name = "v";
        } else {
            v_name = "vg";
        }

        for (int i = 0; i < length(mean_vi); i++){
            model_pars.push_back("v" + std::to_string(k+1) + "." + std::to_string(i+1));
            mean_vi_save.push_back(mean_vi[i]);
        }
        sample_vi_save.conservativeResize(NoChange, sample_vi_save.cols()+sample_vi.cols());
        sample_vi_save.block(0,sample_vi_save.cols()-sample_vi.cols(),iter, sample_vi.cols()) = sample_vi;

        cop_vec_new = copula_type; // no need
        latent_cop_vec_new = latent_copula_type; // no need
    } else {            // Infer v and theta
        if (structfactor == 1 || structfactor == 4) {
            v_name = "v";
            theta_name = "theta";
        } else {
            v_name = "vg";
            v_name = "thetag";
        }

        for (int i = 0; i < t_max; i++) {
            model_pars.push_back("v" + std::to_string(k+1) + "." + std::to_string(i+1));
            mean_vi_save.push_back(mean_vi[i]);
        }
        for (int i = t_max; i < length(mean_vi); i++) {
            model_pars.push_back("theta" + std::to_string(k+1) + "." + std::to_string(i+1-t_max));
            mean_vi_save.push_back(mean_vi[i]);
        }
        sample_vi_save.conservativeResize(NoChange, sample_vi_save.cols()+sample_vi.cols());
        sample_vi_save.block(0,sample_vi_save.cols()-sample_vi.cols(),iter, sample_vi.cols()) = sample_vi;

        if (copselect) {
            copula_type = cop_vec_new;
            latent_copula_type = latent_cop_vec_new;
            // if (structfactor == 2){
            //     latent_copula_type = latent_cop_vec_new;
            // }
            // if (structfactor == 3){
            //     latent_copula_type = latent_cop_vec_new;
            // }

        } else {
            cop_vec_new = copula_type;  // no need
            latent_cop_vec_new = latent_copula_type;  // no need
        }


    }
}


void save_hmc(  std::vector<std::vector<std::string> > parameter_names,
                std::vector<std::vector<double> > parameter_values,
                std::vector<string>& model_pars,
                matrix_d& sample_hmc,
                vector_d& mean_hmc)
{
    int num_sample_var = parameter_names[0].size();
    int iteration = parameter_values.size();

    model_pars.resize(0);
    sample_hmc.resize(iteration, num_sample_var);

    for ( int i = 0; i < num_sample_var; i++) {
        model_pars.push_back( parameter_names[0][i] );
    }
    int count = 0;
    for ( std::vector<double>& param_vector : parameter_values) {
        sample_hmc.row(count) = VectorXd::Map(&param_vector[0], num_sample_var);
        count++;
    }

    mean_hmc = sample_hmc.colwise().mean();
}

void gidtomatrix(int n_max,int n_group,
                 std::vector<int>& gid,
                 std::vector<std::vector<int> >& g_mat,
                 std::vector<int>& g_count){
    for( int i = 0; i < n_max; i++){
        gid[i]--;
        g_mat[gid[i]][g_count[gid[i]]] = i;
        g_count[gid[i]]++;
    }
    for( int i = 0; i < n_group; i++)
    {
        Rcpp::Rcout << " g_count " << i << ": " <<  g_count[i] << std::endl;
    }
}

}
#endif // MISCELLANEOUS_HPP

