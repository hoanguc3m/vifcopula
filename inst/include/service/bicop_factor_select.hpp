// Code generated by Stan version 2.14

#ifndef VIFCOPULA_SERVICE_BICOP_FACTOR_SELECT_HPP
#define VIFCOPULA_SERVICE_BICOP_FACTOR_SELECT_HPP

#include <bicopula_stanc.hpp>
#include <stan/services/optimize/bfgs.hpp>
#include <stan/optimization/bfgs.hpp>


namespace vifcopula{

typedef boost::ecuyer1988 rng_t;
typedef vifcopula::bicopula bicopula;
typedef stan::optimization::BFGSLineSearch<bicopula,stan::optimization::BFGSUpdate_HInv<> > Optimizer_BFGS;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;


double bicop_factor_select(std::vector<double>& u,
                    matrix_d& v2g,
                    int t_max,
                    std::vector<double>& params_out,
                    int& gid_out,
                    rng_t& base_rng){

    std::vector<double> params_r(2);
    params_r[0] = 1;
    params_r[1] = 0;
    std::vector<int> params_i(0);
    bool save_iterations = false;
    int refresh = 0;
    int return_code;
    int return_cop = 0;
    double lpmax = std::numeric_limits<double>::min();
    double BICmin = std::numeric_limits<double>::max();
    int num_factor = v2g.cols();

    // v2g includes v0 at col 0
    for (int j = 1; j < num_factor; j++){
        //v2g_temp = v2g.col(j); get the j_th latent
        vector<double> v2g_temp(t_max);
        VectorXd::Map(&v2g_temp[0], t_max) = v2g.col(j);

        bicopula biuv(1,u,v2g_temp,t_max,base_rng);

        if (biuv.check_Ind()){
            //biuv.set_copula_type(0);
            return return_cop;

        } else {
            const int cop_seq_size = 15;                     // Change the number
            int cop_seq[cop_seq_size] = {1, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36, 2};    // Choose among copula type
            double log_cop[cop_seq_size] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            double AIC[cop_seq_size] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            double BIC[cop_seq_size] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

            int imax=0;
            for (int i = 0; i < cop_seq_size; i++) {
                biuv.set_copula_type(cop_seq[i]);
                std::stringstream out;
                Optimizer_BFGS bfgs(biuv, params_r, params_i, &out);
                double lp = 0;
                int ret = 0;
                while (ret == 0) {
                    ret = bfgs.step();
                }
                lp = bfgs.logp();
                log_cop[i] = lp;


                if (cop_seq[i] == 2){
                    AIC[i] = -2 * lp + 2 * 2;
                    BIC[i] = -2 * lp + log(t_max) * 2;
                } else {
                    AIC[i] = -2 * lp + 2 * 1;
                    BIC[i] = -2 * lp + log(t_max) * 1;
                }

                if (BIC[i] < BICmin){
                    lpmax = lp;
                    BICmin = BIC[i];
                    imax = i;
                    gid_out = j-1;
                    // std::vector<double> get_param;
                    bfgs.params_r(params_r);
                    biuv.write_array(base_rng, params_r, params_i, params_out);
                    return_cop = cop_seq[imax];
                }

            };



        }

    }

    params_out.resize(2); // For theta, theta2 even we dont use it.

    std::cout << " Select cop " << return_cop << " Lp " << lpmax << " "
              << " BIC " << BICmin << " " << params_out[0] << " "
              << params_out[1] << " gid " << gid_out+1 << std::endl ;
    return return_cop;
}   // end func

} // namespace
#endif // VIFCOPULA_SERVICE_BICOP_FACTOR_SELECT_HPP
