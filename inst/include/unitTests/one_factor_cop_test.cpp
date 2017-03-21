#include <gtest/gtest.h>
#include <one_factor_cop.cpp>
#include <iostream>
#include <stan/math.hpp>
#include <boost/exception/all.hpp>
#include <boost/random/additive_combine.hpp>
#include <boost/random/linear_congruential.hpp>

#include <stan/variational/advi.hpp>
#include <stan/interface_callbacks/writer/stream_writer.hpp>

//typedef vifcopula::one_factor_cop Model_cp;

using namespace Eigen;
using namespace stan::math;
//using namespace vifcopula;

typedef Eigen::Matrix<int,Eigen::Dynamic,1> vector_int;
typedef Eigen::Matrix<int,1,Eigen::Dynamic> row_vector_int;
typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> matrix_int;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;
typedef boost::ecuyer1988 rng_t;

//
//TEST(vifobject,  fcopula) {
//
//    // Initiate model
//
//    // Set seed
//    int seed  = 10;
//    rng_t base_rng(seed);
//
//        int t_max  = 10;
//        int n_max  = 3;
//        int k_max  = 1;
//        matrix_d u(t_max,n_max);
//        for (int i = 0; i < t_max; i++ ){
//            for (int j = 0; j < n_max; j++ ){
//                    u(i,j) = uniform_rng(0.0,1.0,base_rng);
//
//            }
//        }
//        std::cout << "u value" << " " << u << std::endl;
//
//        vector_int gid(n_max);
//        gid << 1,1,1;
//
//        vector_d v(t_max);
//        for (int i = 0; i < t_max; i++ ){
//                v(i) = uniform_rng(0.0,1.0,base_rng);
//        }
//
//
//        vector_int copula_type(n_max);
//        copula_type << 0,1,3;
//        k_max = 0;
//
//    // Create object test
//    Model_cp copula_l1(u,gid,copula_type,t_max, n_max, k_max, base_rng);
//        std::cout << "v value" << " " << v << std::endl;
//
//    // Dummy input
//    Eigen::VectorXd cont_params = Eigen::VectorXd::Zero(copula_l1.num_params_r());
//    stan::variational::normal_meanfield cont_params_mf(cont_params);
//
//    //stan::interface_callbacks::writer::stream_writer writer;
//
//    stan::variational::advi<Model_cp, stan::variational::normal_meanfield, rng_t> test_advi(copula_l1,
//                                                                                            cont_params,
//                                                                                            base_rng,
//                                                                                            10,
//                                                                                            300,
//                                                                                            100,
//                                                                                            1000);
//
//    std::cout << " Number of Parameters:" << copula_l1.num_params_r() << std::endl;
//    //std::cout << " calc ELBO :" << test_advi.calc_ELBO(cont_params_mf,writer) << std::endl;
//
//}
