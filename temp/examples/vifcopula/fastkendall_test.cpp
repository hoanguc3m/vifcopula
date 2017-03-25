#include <gtest/gtest.h>
#include <BiCopIndTest.hpp>
#include <boost/random/additive_combine.hpp> // L'Ecuyer RNG
#include <stan/math.hpp>

typedef boost::ecuyer1988 rng_t;

typedef Eigen::Matrix<int,Eigen::Dynamic,1> vector_int;
typedef Eigen::Matrix<int,1,Eigen::Dynamic> row_vector_int;
typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> matrix_int;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;

using stan::math::uniform_rng;
using std::vector;


TEST(advi_test, kendall_test) {
          // Set seed
    int seed  = 10;
    rng_t base_rng(seed);

        int t_max  = 10;
        int n_max  = 1;
        int k_max  = 1;
        vector<double> u(t_max);
        for (int i = 0; i < t_max; i++ ){
            u[i] = uniform_rng(0.0,1.0,base_rng);
        }

        vector<int> gid = {1,1,1,1,1,1,1};


        vector<double> v(t_max);
        for (int i = 0; i < t_max; i++ ){
                v[i] = uniform_rng(0.0,1.0,base_rng);
        }

        vector<int> copula_type = {0,1,2,3,4,5,6};
        k_max = 0;

    double p_val = BiCopIndTest(u,v);
    std::cout << " p.val = " << p_val << std::endl;

}
