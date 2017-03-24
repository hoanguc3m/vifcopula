#include <one_factor_cop.hpp>
#include <advi_mod.hpp>
#include <stan/interface_callbacks/writer/stream_writer.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>
#include <string>
#include <iostream>
#include <boost/random/additive_combine.hpp> // L'Ecuyer RNG
#include <stan/model/log_prob_propto.hpp>
#include <stan/model/grad_hess_log_prob.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using stan::math::uniform_rng;
using std::vector;


typedef boost::ecuyer1988 rng_t;
typedef vifcopula::one_factor_cop Model_cp;
typedef Eigen::Matrix<int,Eigen::Dynamic,1> vector_int;
typedef Eigen::Matrix<int,1,Eigen::Dynamic> row_vector_int;
typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> matrix_int;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;



TEST(advi_test, one_factor_cop_constraint_meanfield) {
  // Create mock data_var_context
  std::fstream data_stream("hier_logistic.data.R",
                           std::fstream::in);
  stan::io::dump data_var_context(data_stream);
  data_stream.close();

  std::stringstream output;
  output.clear();

      // Set seed
    int seed  = 10;
    rng_t base_rng(seed);

        int t_max  = 10;
        int n_max  = 7;
        int k_max  = 1;
        matrix_d u(t_max,n_max);
        for (int i = 0; i < t_max; i++ ){
            for (int j = 0; j < n_max; j++ ){
                    u(i,j) = uniform_rng(0.0,1.0,base_rng);

            }
        }
        std::cout << "u value" << " " << u << std::endl;

        vector<int> gid = {1,1,1,1,1,1,1};


        vector_d v(t_max);
        for (int i = 0; i < t_max; i++ ){
                v(i) = uniform_rng(0.0,1.0,base_rng);
        }


        vector<int> copula_type = {0,1,2,3,4,5,6};
        k_max = 0;

  std::vector<double> params_r(t_max+n_max);
  std::vector<int> params_i(0);
  std::vector<double> gradient;
  std::vector<double> hessian;

  // Instantiate model
  Model_cp my_model(u,gid,copula_type,t_max, n_max, k_max, base_rng);

  // RNG
  //rng_t base_rng(0);

  // Dummy input
    std::cout << " copula LL :" << my_model.num_params_r() << std::endl;

  Eigen::VectorXd cont_params = Eigen::VectorXd::Zero(my_model.num_params_r());
    for (int i = 0; i < t_max; i++ ){
                cont_params(i) = v(i);
                params_r[i] = v(i);
        }
  cont_params(t_max) = 0.5;cont_params(t_max+1) = 0.5;cont_params(t_max+2) = 5;
  cont_params(t_max+3) = 0.5;cont_params(t_max+4) = 1.5;cont_params(t_max+5) = 0.5;
  cont_params(t_max+6) = 1.5;


  stan::variational::normal_meanfield cont_params1(t_max + n_max);


  // ADVI
  stan::variational::advi_mod<Model_cp, stan::variational::normal_meanfield, rng_t> test_advi(my_model,
                                                     cont_params,
                                                     base_rng,
                                                     10,
                                                     100,
                                                     10,
                                                     2);
  std::stringstream out_message_writer;
  stan::interface_callbacks::writer::stream_writer message_writer(std::cout);

  std::stringstream out_parameter_writer;
  stan::interface_callbacks::writer::stream_writer parameter_writer(out_parameter_writer);

  std::stringstream out_diagnostic_writer;
  stan::interface_callbacks::writer::stream_writer diagnostic_writer(std::cout);


  std::cout << " copula LL :" << stan::model::log_prob_propto<true, stan_model>(my_model,cont_params,0) << std::endl;

//  params_r[t_max] = 0;params_r[t_max+1] = 0.5;params_r[t_max+2] = 0.5;
//  params_r[t_max+3] = 0.5;params_r[t_max+4] = 1.5;params_r[t_max+5] = 0.5;
//  params_r[t_max+6] = 1.5;
//
//  stan::model::grad_hess_log_prob<false, true, stan_model>(my_model, params_r, params_i, gradient, hessian, 0);
//
//  std::cout << " gradient :" << gradient[0] << " " << gradient[1] << std::endl;
//  std::cout << " hessian :" << hessian[0] << " " << hessian[1] << std::endl;

  test_advi.run(0.01, true, 50, 1, 2e4,
                message_writer, parameter_writer, diagnostic_writer);

   std::cout << " out_stream " << out_parameter_writer.str() << std::endl;


    int iter = 2;



    int max_param = my_model.num_params_r();
    matrix_d sample_iv(iter,max_param);
    vector_d mean_iv(max_param);

    std::string token;
    std::getline(out_parameter_writer, token);
    std::getline(out_parameter_writer, token);

    int j;
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

        for (int i = 0; i < iter;i++){
            std::getline(out_parameter_writer, token, ',');
                for (int j = 0; j < max_param-1;j++){
                    std::getline(out_parameter_writer, token, ',');
                    boost::trim(token);
                    sample_iv(i,j) = boost::lexical_cast<double>(token);
                }
            std::getline(out_parameter_writer, token);
                boost::trim(token);
                sample_iv(i,max_param-1) = boost::lexical_cast<double>(token);
        }

    } catch (...){
        std::cout << " j " << j << std::endl;
        std::cout << " token " << token << std::endl;
    }







   std::cout << " sample_iv " << sample_iv << std::endl;

}


