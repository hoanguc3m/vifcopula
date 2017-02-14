#ifndef VIFCOPULA_BICOPDIST_H
#define VIFCOPULA_BICOPDIST_H

#include <stan/math.hpp>

namespace vifcopula {

    template <bool propto, typename T_u, typename T_v>
    typename stan::return_type<T_u, T_v>::type
        bicop_independence_log(const T_u& u, const T_v& v);

    template <typename T_u, typename T_v>
    inline
        typename stan::return_type<T_u, T_v>::type
        bicop_independence_log(const T_u& u, const T_v& v);

    template <bool propto, typename T_u, typename T_v, typename T_rho>
    typename stan::return_type<T_u, T_v, T_rho>::type
        bicop_normal_log(const T_u& u, const T_v& v, const T_rho& rho);

    template <typename T_u, typename T_v, typename T_rho>
    inline
        typename stan::return_type<T_u, T_v, T_rho>::type
        bicop_normal_log(const T_u& u, const T_v& v, const T_rho& rho);



}
#endif // VIFCOPULA_BICOPDIST_H
