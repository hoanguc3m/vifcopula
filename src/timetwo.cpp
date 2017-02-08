#include <Rcpp.h>
#include <RcppEigen.h>
#include <stan/math.hpp>

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(StanHeaders)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
int sum_example() {

    using std::pow;
    std::vector<double> y = {1.3,2,3};
    std::vector<double> x = {2,4,6};

    stan::math::var beta = 0.5, sigma = 1.2;
    std::vector<stan::math::var> acc;

    for (int n = 0; n < y.size(); ++n)
        acc.push_back(stan::math::normal_log(y[n], x[n] * beta, sigma));
    stan::math::var lp = stan::math::sum(acc);

    lp.grad();
    std::cout << " d.f / d.beta = " << beta.adj()
              << " d.f / d.sigma = " << sigma.adj() << std::endl;

    return 0;
}


