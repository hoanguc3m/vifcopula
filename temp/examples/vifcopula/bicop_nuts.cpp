#include <bicopula_stanc.hpp>

#include <stan/callbacks/stream_writer.hpp>
#include <stan/callbacks/stream_logger.hpp>
#include <stan/callbacks/interrupt.hpp>

#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <iostream>
#include <boost/random/additive_combine.hpp> // L'Ecuyer RNG
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <stan/services/sample/hmc_nuts_diag_e_adapt.hpp>
#include <stan/io/empty_var_context.hpp>

#include <test/unit/services/instrumented_callbacks.hpp>


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


TEST(nuts_test, bicopula)
{

    std::stringstream output;
    output.clear();

    // Set seed
    int seed  = 10;
    rng_t base_rng(seed);

    int t_max  = 5;
    int k_max  = 1;
    vector<double> u = {0.1,0.2,0.3,0.4,0.5};
//    vector<double> u(t_max);
//    for (int i = 0; i < t_max; i++ )
//    {
//        u[i] = uniform_rng(0.0,1.0,base_rng);
//    }

    vector<int> gid = {1,1,1,1,1,1,1};
    vector<double> v = {0.25,0.2,0.35,0.2,0.3};

//    vector<double> v(t_max);
//    for (int i = 0; i < t_max; i++ )
//    {
//        //v[i] = (u[i] + uniform_rng(0.0,1.0,base_rng))/2;
//        v[i] = uniform_rng(0.0,1.0,base_rng);
//
//    }
    int copula_type = 2;
    k_max = 0;

    std::vector<double> params_r(2);
    std::vector<int> params_i(0);
    std::vector<double> gradient;
    std::vector<double> hessian;

    // Instantiate model
    Model_cp my_model(copula_type,u,v,t_max, base_rng);

    // RNG
    //rng_t base_rng(0);

    // Dummy input
    std::cout << " copula params :" << my_model.num_params_r() << std::endl;

    Eigen::VectorXd cont_params = Eigen::VectorXd::Zero(my_model.num_params_r());
    cont_params(0) = 1;
    params_r[0] = 1;
    params_r[1] = 0;

    // stan::callbacks::stream_logger logger(debug, info, warn, error, fatal);
    stan::callbacks::stream_logger logger(std::cout,std::cout,std::cout,std::cout,std::cout);
    //stan::test::unit::instrumented_logger logger;


//    std::stringstream parameter_stream_;
//    std::stringstream diagnostic_stream_;
//    stan::callbacks::stream_writer init_writer(diagnostic_stream_), sample_writer(parameter_stream_), diagnostic_writer(std::cout);

    stan::test::unit::instrumented_writer init_writer, sample_writer, diagnostic_writer;


    unsigned int random_seed = 0;
    unsigned int chain = 1;
    double init_radius = 0;
    int num_warmup = 1000;
    int num_samples = 2000;
    int num_thin = 1;
    bool save_warmup = true;
    int refresh = (num_warmup + num_samples) /10;
    double stepsize = 1;
    double stepsize_jitter = 0;
    int max_depth = 10;
    double delta = 0.80000000000000004;
    double gamma = 0.050000000000000003;
    double kappa = .75;
    double t0 = 10;
    unsigned int init_buffer = 75;
    unsigned int term_buffer = 50;
    unsigned int window = 25;


    //stan::callbacks::interrupt interrupt;
    stan::test::unit::instrumented_interrupt interrupt;

//    int return_code = stan::services::sample::hmc_nuts_diag_e_adapt(
//                          my_model, params_r, random_seed, chain, init_radius,
//                          num_warmup, num_samples, num_thin, save_warmup, refresh,
//                          stepsize, stepsize_jitter, max_depth, delta, gamma, kappa, t0,
//                          init_buffer, term_buffer, window,
//                          interrupt, logger, init,
//                          parameter, diagnostic);

        boost::ecuyer1988 rng = stan::services::util::create_rng(random_seed, chain);

        std::vector<int> disc_vector;
        std::vector<double> cont_vector
          = params_r;

        Eigen::VectorXd inv_metric;

        stan::io::dump dmp = stan::services::util::create_unit_e_diag_inv_metric(my_model.num_params_r());
        stan::io::var_context& unit_e_metric = dmp;

        try {
          inv_metric =
            stan::services::util::read_diag_inv_metric(unit_e_metric, my_model.num_params_r(),
                                        logger);
          stan::services::util::validate_diag_inv_metric(inv_metric, logger);
        } catch (const std::domain_error& e) {
        }

        stan::mcmc::adapt_diag_e_nuts<Model_cp, boost::ecuyer1988>
          sampler(my_model, rng);

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

        stan::services::util::run_adaptive_sampler(sampler, my_model, cont_vector, num_warmup,
                                   num_samples, num_thin, refresh, save_warmup,
                                   rng, interrupt, logger,
                                   sample_writer, diagnostic_writer);


        std::vector<std::vector<std::string> > parameter_names;
        parameter_names = sample_writer.vector_string_values();
        std::vector<std::vector<double> > parameter_values;
        parameter_values = sample_writer.vector_double_values();

        // Expectations of parameter parameter names.
        ASSERT_EQ(9, parameter_names[0].size());
        EXPECT_EQ("lp__", parameter_names[0][0]);
        EXPECT_EQ("accept_stat__", parameter_names[0][1]);
        EXPECT_EQ("stepsize__", parameter_names[0][2]);
        EXPECT_EQ("treedepth__", parameter_names[0][3]);
        EXPECT_EQ("n_leapfrog__", parameter_names[0][4]);
        EXPECT_EQ("divergent__", parameter_names[0][5]);
        EXPECT_EQ("energy__", parameter_names[0][6]);
        EXPECT_EQ("theta", parameter_names[0][7]);
        EXPECT_EQ("theta2", parameter_names[0][8]);

        // Expect one name per parameter value.
        EXPECT_EQ(parameter_names[0].size(), parameter_values[0].size());
        EXPECT_EQ((num_warmup+num_samples)/num_thin, parameter_values.size());


        EXPECT_EQ(num_warmup + num_samples, interrupt.call_count());


}



