#ifndef VIFCOPULA_LOGBIFCOP_H
#define VIFCOPULA_LOGBIFCOP_H

namespace vifcopula {
    typedef Eigen::Matrix<int,Eigen::Dynamic,1> vector_int;
    typedef Eigen::Matrix<int,1,Eigen::Dynamic> row_vector_int;
    typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> matrix_int;
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
    typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;

    double logBifcop(const matrix_d& u_, const matrix_d& v_, const matrix_d& par_,
                    const matrix_int& copula_type_, const int& t_max_,const int& i,
                    const int& k_);
}
#endif // VIFCOPULA_LOGBIFCOP_H
