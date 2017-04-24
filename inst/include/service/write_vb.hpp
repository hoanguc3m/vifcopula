#ifndef VIFCOPULA_SERVICE_WRITEVB_CPP
#define VIFCOPULA_SERVICE_WRITEVB_CPP

#include <stan/math.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;
typedef Eigen::Matrix<int,Eigen::Dynamic,1> vector_int;
typedef Eigen::Matrix<int,1,Eigen::Dynamic> row_vector_int;
typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> matrix_int;

/* PRINT_ELEMENTS()
 * - prints optional C-string optcstr followed by
 * - all elements of the collection coll
 * - separated by spaces
 */
template <class T>
inline void PRINT_ELEMENTS (const T& coll, const char* optcstr="")
{
    typename T::const_iterator pos;

    std::cout << optcstr;
    for (pos=coll.begin(); pos!=coll.end(); ++pos) {
        std::cout << *pos << ' ';
    }
    std::cout << std::endl;
}


void write_vb (std::stringstream& out_parameter_writer, vector_d& mean_iv, matrix_d& sample_iv){

    std::string token;
    // std::getline(out_parameter_writer, token);
    // std::getline(out_parameter_writer, token);
    int iter = sample_iv.rows();
    int max_param = sample_iv.cols();
    int i,j;
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

        for (i = 0; i < iter;i++){
            std::getline(out_parameter_writer, token, ',');
            for (j = 0; j < max_param-1;j++){
                std::getline(out_parameter_writer, token, ',');
                boost::trim(token);
                sample_iv(i,j) = boost::lexical_cast<double>(token);
            }
            std::getline(out_parameter_writer, token);
            boost::trim(token);
            sample_iv(i,max_param-1) = boost::lexical_cast<double>(token);
        }

    } catch (...){
        std::cout << "Error at i = " << i << " j = " << j << std::endl;
        std::cout << " token " << token << std::endl;
    }

}


#endif // VIFCOPULA_SERVICE_WRITEVB_CPP
