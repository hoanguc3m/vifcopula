#ifndef VIFCOPULA_HFUNC_HPP
#define VIFCOPULA_HFUNC_HPP

#include <transform/hfunc_independence.hpp>
#include <transform/hfunc_normal.hpp>
#include <transform/hfunc_student.hpp>
#include <transform/hfunc_clayton.hpp>
#include <transform/hfunc_gumbel.hpp>
#include <transform/hfunc_frank.hpp>
#include <transform/hfunc_joe.hpp>

namespace vifcopula
{

using namespace stan::math;
using namespace vifcopula;

typedef Eigen::Matrix<int,Eigen::Dynamic,1> vector_int;
typedef Eigen::Matrix<int,1,Eigen::Dynamic> row_vector_int;
typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> matrix_int;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;


void hfunc_trans(matrix_d& u, const vector_d& mean_iv, const std::vector<int>& cop_vec){
    int t_max = u.rows();
    int n_max = u.cols();
    int theta_pos = t_max;
    double theta1;
    double theta2;

    vector_d v_temp = mean_iv.head(t_max);
    vector_d u_temp(t_max);
    vector_d u_trans(t_max);

    for (int i = 0; i < n_max; i++){
        u_temp = u.col(i);
        switch ( cop_vec[i] ) {
        case 0:
            // Independence copula
            u_trans = hfunc_independence(u_temp,v_temp);
            break;
        case 1:
            // Gaussian copula
            theta1 = mean_iv[theta_pos]; theta_pos++;
            u_trans = hfunc_normal(u_temp,v_temp,theta1);
            break;
        case 2:
            // Student copula
            theta1 = mean_iv[theta_pos]; theta_pos++;
            theta2 = mean_iv[theta_pos]; theta_pos++;
            u_trans = hfunc_student(u_temp,v_temp,theta1,theta2);
            break;
        case 3:
            // Clayon copula
            theta1 = mean_iv[theta_pos]; theta_pos++;
            u_trans = hfunc_clayton(u_temp,v_temp,theta1);
            break;
        case 4:
            // Gumbel copula
            theta1 = mean_iv[theta_pos]; theta_pos++;
            u_trans = hfunc_gumbel(u_temp,v_temp,theta1);
            break;
        case 5:
            // Frank copula
            theta1 = mean_iv[theta_pos]; theta_pos++;
            u_trans = hfunc_frank(u_temp,v_temp,theta1);
            break;
        case 6:
            // Joe copula
            theta1 = mean_iv[theta_pos]; theta_pos++;
            u_trans = hfunc_joe(u_temp,v_temp,theta1);
            break;
        default:
            // Code to execute if <variable> does not equal the value following any of the cases
            // Send a break message.
            break;
        }
        u.col(i) = u_trans;

    }   // end for h_transformation
}       // end function

}       // end name_sapce
#endif // VIFCOPULA_HFUNC_HPP


