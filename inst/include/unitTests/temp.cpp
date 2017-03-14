#ifndef VIFCOPULA_DISTRIBUTION_TEMP_INDEPENDENCE_LOG_TEST_CPP
#define VIFCOPULA_DISTRIBUTION_TEMP_INDEPENDENCE_LOG_TEST_CPP

#include <stan/math.hpp>
#include <dist/bicop_normal_log.cpp>

#include <ostream>

void tempfunction(void) {
    using stan::math::var;
    for (double v_val = -.5; v_val = 0; v_val = 0.5) {
        var v(v_val);
        var lp1(0.0);
        lp1 += vifcopula::bicop_normal_log<true>(0.1, v, 0.5);
        double lp1val = lp1.val();
        stan::math::grad(lp1.vi_);
        double lp1adj = lp1.adj();
        std::cout << lp1adj << " " << lp1adj << std::endl;

    }
}

#endif // VIFCOPULA_DISTRIBUTION_TEMP_INDEPENDENCE_LOG_TEST_CPP
