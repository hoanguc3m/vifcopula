#include <bicopula.hpp>
#include <advi_mod.hpp>
#include <stan/services/optimize/do_bfgs_optimize.hpp>
#include <stan/optimization/bfgs.hpp>

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
typedef vifcopula::bicopula Model_cp;
typedef Eigen::Matrix<int,Eigen::Dynamic,1> vector_int;
typedef Eigen::Matrix<int,1,Eigen::Dynamic> row_vector_int;
typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> matrix_int;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;

struct mock_callback
{
    int n;
    mock_callback() : n(0) { }

    void operator()()
    {
        n++;
    }
};

TEST(advi_test, bicopula)
{
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

    int t_max  = 100;
    int k_max  = 1;
    vector<double> u(t_max);
    for (int i = 0; i < t_max; i++ )
    {
        u[i] = uniform_rng(0.0,1.0,base_rng);
    }

    vector<int> gid = {1,1,1,1,1,1,1};

    vector<double> v(t_max);
    for (int i = 0; i < t_max; i++ )
    {
        v[i] = uniform_rng(0.0,1.0,base_rng);
    }
    int copula_type = 6;
    k_max = 0;

    std::vector<double> params_r(1);
    std::vector<int> params_i(0);
    std::vector<double> gradient;
    std::vector<double> hessian;

    // Instantiate model
    Model_cp my_model(copula_type,u,v,t_max, base_rng);

    // RNG
    //rng_t base_rng(0);

    // Dummy input
    std::cout << " copula LL :" << my_model.num_params_r() << std::endl;

    Eigen::VectorXd cont_params = Eigen::VectorXd::Zero(my_model.num_params_r());
    cont_params(0) = 1;
    params_r[0] = 1;

    /*
          // ADVI
          stan::variational::advi_mod<Model_cp, stan::variational::normal_meanfield, rng_t> bicop(my_model,
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
          stan::interface_callbacks::writer::stream_writer diagnostic_writer(out_diagnostic_writer);


          bicop.run(0.01, true, 50, 0.1, 2e4,
                        message_writer, parameter_writer, diagnostic_writer);

           std::cout << " out_stream " << out_parameter_writer.str() << std::endl;

    */

    // do_bfgs_optimize__bfgs
    typedef stan::optimization::BFGSLineSearch<Model_cp,stan::optimization::BFGSUpdate_HInv<> > Optimizer_BFGS;
    std::stringstream out;
    Optimizer_BFGS bfgs(my_model, params_r, params_i, &out);
    std::cout << " out " << out.str() << std::endl;

    double lp = 0;
    bool save_iterations = true;
    int refresh = 0;
    int return_code;

    mock_callback callback;

    stan::interface_callbacks::writer::stream_writer writer(out);
    std::stringstream info_ss;
    stan::interface_callbacks::writer::stream_writer info(info_ss);
    return_code = stan::services::optimize::do_bfgs_optimize(my_model, bfgs, base_rng,
                  lp, params_r, params_i,
                  writer, info,
                  save_iterations, refresh,
                  callback);

    std::cout << " out_stream " << info_ss.str() << std::endl;
    std::cout << " return_code " << return_code << std::endl;
    std::cout << " callback.n " << callback.n << std::endl;
    std::cout << " lp " << lp << std::endl;
    std::cout << " writer " << out.str() << std::endl;
    std::cout << " Independent test " << my_model.check_Ind() << std::endl;
}



