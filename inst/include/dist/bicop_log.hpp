#ifndef VIFCOPULA_DISTRIBUTION_BICOP_LOG_HPP
#define VIFCOPULA_DISTRIBUTION_BICOP_LOG_HPP

#include <dist/bicop_independence_log.hpp>
#include <dist/bicop_normal_log.hpp>
#include <dist/bicop_student_log.hpp>
#include <dist/bicop_clayton_log.hpp>
#include <dist/bicop_gumbel_log.hpp>
#include <dist/bicop_frank_log.hpp>
#include <dist/bicop_joe_log.hpp>

#include <dist/bicop_survival_clayton_log.hpp>
#include <dist/bicop_survival_gumbel_log.hpp>
#include <dist/bicop_survival_joe_log.hpp>

#include <dist/bicop_r90_clayton_log.hpp>
#include <dist/bicop_r90_gumbel_log.hpp>
#include <dist/bicop_r90_joe_log.hpp>

#include <dist/bicop_r270_clayton_log.hpp>
#include <dist/bicop_r270_gumbel_log.hpp>
#include <dist/bicop_r270_joe_log.hpp>

namespace vifcopula {

using namespace stan::math;
using namespace stan;

template <bool propto__, bool jacobian__,
                typename Tu__, typename Tv__,typename T__ >

void bicop_log_add(int i,
    const std::vector<int>& copula_type,
    Eigen::Matrix<Tu__,Eigen::Dynamic,1>& u,
    Eigen::Matrix<Tv__,Eigen::Dynamic,1>&  v,
    Eigen::Matrix<T__,Eigen::Dynamic,1>& theta,
    Eigen::Matrix<T__,Eigen::Dynamic,1>&  theta2,
    T__& lp__,
    stan::math::accumulator<T__>& lp_accum__,
    stan::io::reader<T__>& in__
    ) {
    static const char* function("vifcopula::bicop_log_add");

    int ibase = i+1;
    switch ( copula_type[i] ) {
    case 0:
        // Independence copula
        //theta[i] = in__.scalar_constrain();
        theta[i] = 0;
        lp_accum__.add(bicop_independence_log<propto__>(u,v));
        break;
    case 1:
        // Gaussian copula
        // if (i == 0){
        //     if (jacobian__)
        //         theta[i] = in__.scalar_lub_constrain(0,1,lp__);
        //     else
        //         theta[i] = in__.scalar_lub_constrain(0,1);
        // } else {
        //     if (jacobian__)
        //         theta[i] = in__.scalar_lub_constrain(-0.9,1,lp__);
        //     else
        //         theta[i] = in__.scalar_lub_constrain(-0.9,1);
        // }

        if (jacobian__)

            theta[i] = in__.scalar_lub_constrain(0,1,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(0,1);

        lp_accum__.add(bicop_normal_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        break;
    case 2:
        // Student copula

        // if (i == 0){
        //     if (jacobian__)
        //         theta[i] = in__.scalar_lub_constrain(0,1,lp__);
        //     else
        //         theta[i] = in__.scalar_lub_constrain(0,1);
        // } else {
        //     if (jacobian__)
        //         theta[i] = in__.scalar_lub_constrain(-0.9,1,lp__);
        //     else
        //         theta[i] = in__.scalar_lub_constrain(-0.9,1);
        // }

        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(0,1,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(0,1);


        if (jacobian__)
            theta2[i] = in__.scalar_lub_constrain(2.1,40,lp__);
        else
            theta2[i] = in__.scalar_lub_constrain(2.1,40);

        lp_accum__.add(bicop_student_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1),
            get_base1(theta2,ibase,"theta",1)));

        // lp_accum__.add(bicop_student_log<propto__>(u,
        //     v,
        //     get_base1(theta,ibase,"theta",1),
        //     7));
        break;
    case 3:
        // Clayon copula

        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(0.001,30,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(0.001,30);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_clayton_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        break;
    case 4:
        // Gumbel copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(1,15,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(1,15);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 1, Inf)); //Improper priors
        lp_accum__.add(bicop_gumbel_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        break;
    case 5:
        // Frank copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(0.1,100,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(0.1,100);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_frank_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));

        break;
    case 6:
        // Joe copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(1,30,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(1,30);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_joe_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));

        break;

    case 13:
        // survival Clayon copula

        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(0.001,30,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(0.001,30);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_survival_clayton_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        break;
    case 14:
        // survival Gumbel copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(1,15,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(1,15);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 1, Inf)); //Improper priors
        lp_accum__.add(bicop_survival_gumbel_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        break;
    case 16:
        // survival Joe copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(1,30,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(1,30);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_survival_joe_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));

        break;

    case 21:

        if (jacobian__)

            theta[i] = in__.scalar_lub_constrain(-1,0,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-1,0);

        lp_accum__.add(bicop_normal_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        break;
    case 22:
        // Student copula

        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-1,0,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-1,0);

        if (jacobian__)
            theta2[i] = in__.scalar_lub_constrain(2.1,40,lp__);
        else
            theta2[i] = in__.scalar_lub_constrain(2.1,40);

        lp_accum__.add(bicop_student_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1),
            get_base1(theta2,ibase,"theta",1)));

        break;

    case 23:
        // rotated 90 degree Clayon copula

        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-30,-0.001,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-30,-0.001);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_r90_clayton_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        break;
    case 24:
        // rotated 90 degree Gumbel copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-15,-1,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-15,-1);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 1, Inf)); //Improper priors
        lp_accum__.add(bicop_r90_gumbel_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        break;

    case 25:
        // Frank copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-100,-0.1,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-100,-0.1);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_frank_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));

        break;

    case 26:
        // rotated 90 degree Joe copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-30,-1,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-30,-1);
        // std::cout << " theta " << theta[i] << " " << min(v) << " " << max(v) << std::endl;
        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_r90_joe_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));

        break;

    case 33:
        // rotated 270 degree Clayon copula

        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-30,-0.001,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-30,-0.001);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_r270_clayton_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        break;
    case 34:
        // rotated 270 degree Gumbel copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-15,-1,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-15,-1);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 1, Inf)); //Improper priors
        lp_accum__.add(bicop_r270_gumbel_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        break;
    case 36:
        // rotated 270 degree Joe copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-30,-1,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-30,-1);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_r270_joe_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));

        break;

    default:
        // Code to execute if <variable> does not equal the value following any of the cases
        // Send a break message.
        break;
    }



}


}




#endif // VIFCOPULA_DISTRIBUTION_BICOP_LOG_HPP


