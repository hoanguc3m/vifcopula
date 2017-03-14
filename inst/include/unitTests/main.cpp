#include "gtest/gtest.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

//#include <iostream>
//#include <dist/bicop_normal_log.cpp>
//
//int main() {
//    using stan::math::var;
//     std::cout << " Hello world " << std::endl;
//
//
//        double v_val = .5;
//        var v(v_val);
//        var lp1(0.0);
//        lp1 += vifcopula::bicop_normal_log<true>(0.1, v, 0.5);
//        double lp1val = lp1.val();
//        stan::math::grad(lp1.vi_);
//        double lp1adj = lp1.adj();
//        std::cout << lp1adj << " " << lp1adj << std::endl;
//
//
//  return 0;
//}
