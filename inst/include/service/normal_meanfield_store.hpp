#ifndef VIFCOPULA_SERVICE_NORMAL_MEANFIELD_STORE_HPP
#define VIFCOPULA_SERVICE_NORMAL_MEANFIELD_STORE_HPP

#include <stan/math.hpp>

namespace vifcopula
{
    struct normal_meanfield_store{
        /**
         * Mean vector.
         */
        Eigen::VectorXd mu_;

        /**
         * Log standard deviation (log scale) vector.
         */
        Eigen::VectorXd omega_;
    };
}// namespace
#endif // VIFCOPULA_SERVICE_NORMAL_MEANFIELD_STORE_HPP
