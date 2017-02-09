#include <Rcpp.h>
#include <RcppEigen.h>
#include <stan/math.hpp>


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(StanHeaders)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::interfaces(r, cpp)]]

using namespace Rcpp;
using namespace Eigen;
using namespace stan::math;

//' Sum of vector elements.
//'
//' \code{sum_example} returns the sum of all the values present in its arguments.
//'
//' This is a generic function: methods can be defined for it directly
//' or via the \code{\link{Summary}} group generic. For this to work properly,
//' the arguments \code{...} should be unnamed, and dispatch is on the
//' first argument.
//'
//' @param ... Numeric, complex, or logical vectors.
//' @param na.rm A logical scalar. Should missing values (including NaN)
//'   be removed?
//' @return If all inputs are integer and logical, then the output
//'   will be an integer. If integer overflow
//'   \url{http://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
//'   will be NA with a warning. Otherwise it will be a length-one numeric or
//'   complex vector.
//'
//'   Zero-length vectors have sum 0 by definition. See
//'   \url{http://en.wikipedia.org/wiki/Empty_sum} for more details.
//' @examples
//' sum(1:10)
//' sum(1:5, 6:10)
//' sum(F, F, F, T, T)
//'
//' sum(.Machine$integer.max, 1L)
//' sum(.Machine$integer.max, 1)
//'
//' \dontrun{
//' sum("a")
//' }
//' @export
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


