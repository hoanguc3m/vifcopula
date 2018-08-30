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
    using stan::math::square;
    using stan::math::log;

    double log_2 = log(2);
    double log_pi = log(pi());

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

        if (jacobian__)

            theta[i] = in__.scalar_lub_constrain(0,0.99,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(0,0.99);

        lp_accum__.add(bicop_normal_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( log_2 - log_pi - 0.5 * log(1 - square(theta[i])) );
        break;
    case 2:
        // Student copula

        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(0,0.99,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(0,0.99);


        if (jacobian__)
            theta2[i] = in__.scalar_lub_constrain(2.001,30,lp__);
        else
            theta2[i] = in__.scalar_lub_constrain(2.001,30);

        lp_accum__.add(bicop_student_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1),
            get_base1(theta2,ibase,"theta",1)));
        lp_accum__.add( log_2 - log_pi - 0.5 * log(1 - square(theta[i])) );
        lp_accum__.add( gamma_lpdf(theta[i], 1, 0.1)  );
        break;
    case 3:
        // Clayton copula

        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(0.001,20,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(0.001,20);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_clayton_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( log_2 - 2 * log( theta[i] + 2 ) );

        break;
    case 4:
        // Gumbel copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(1,10,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(1,10);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 1, Inf)); //Improper priors
        lp_accum__.add(bicop_gumbel_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( - 2 * log(theta[i]) );

        break;
    case 5:
        // Frank copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(0.1,40,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(0.1,40);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_frank_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( cauchy_lpdf(theta[i], 0, 6) );

        break;
    case 6:
        // Joe copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(1,20,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(1,20);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_joe_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( log_2 - 2 * log( theta[i] + 2 ) );

        break;

    case 13:
        // survival Clayton copula

        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(0.001,20,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(0.001,20);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_survival_clayton_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( log_2 - 2 * log( theta[i] + 2 ) );
        break;
    case 14:
        // survival Gumbel copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(1,10,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(1,10);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 1, Inf)); //Improper priors
        lp_accum__.add(bicop_survival_gumbel_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( - 2 * log(theta[i]) );
        break;
    case 16:
        // survival Joe copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(1,20,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(1,20);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_survival_joe_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( log_2 - 2 * log( theta[i] + 2 ) );
        break;

    case 21:

        if (jacobian__)

            theta[i] = in__.scalar_lub_constrain(-0.99,0,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-0.99,0);

        lp_accum__.add(bicop_normal_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( log_2 - log_pi - 0.5 * log(1 - square(theta[i])) );
        break;
    case 22:
        // Student copula

        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-0.99,0,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-0.99,0);

        if (jacobian__)
            theta2[i] = in__.scalar_lub_constrain(2.001,30,lp__);
        else
            theta2[i] = in__.scalar_lub_constrain(2.001,30);

        lp_accum__.add(bicop_student_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1),
            get_base1(theta2,ibase,"theta",1)));
        lp_accum__.add( log_2 - log_pi - 0.5 * log(1 - square(theta[i])) );
        break;

    case 23:
        // rotated 90 degree Clayton copula

        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-20,-0.001,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-20,-0.001);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_r90_clayton_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( log_2 - 2 * log( -theta[i] + 2 ) );
        break;
    case 24:
        // rotated 90 degree Gumbel copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-10,-1,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-10,-1);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 1, Inf)); //Improper priors
        lp_accum__.add(bicop_r90_gumbel_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( - 2 * log(-theta[i]) );
        break;

    case 25:
        // Frank copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-40,-0.1,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-40,-0.1);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_frank_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( cauchy_lpdf(theta[i], 0, 6) );
        break;

    case 26:
        // rotated 90 degree Joe copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-20,-1,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-20,-1);
        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_r90_joe_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( log_2 - 2 * log( -theta[i] + 2 ) );

        break;

    case 33:
        // rotated 270 degree Clayton copula

        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-20,-0.001,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-20,-0.001);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_r270_clayton_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( log_2 - 2 * log( -theta[i] + 2 ) );
        break;
    case 34:
        // rotated 270 degree Gumbel copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-10,-1,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-10,-1);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 1, Inf)); //Improper priors
        lp_accum__.add(bicop_r270_gumbel_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( - 2 * log(-theta[i]) );
        break;
    case 36:
        // rotated 270 degree Joe copula
        if (jacobian__)
            theta[i] = in__.scalar_lub_constrain(-20,-1,lp__);
        else
            theta[i] = in__.scalar_lub_constrain(-20,-1);

        //lp_accum__.add(uniform_lpdf<propto__>(theta[i], 0, Inf)); //Improper priors
        lp_accum__.add(bicop_r270_joe_log<propto__>(u,
            v,
            get_base1(theta,ibase,"theta",1)));
        lp_accum__.add( log_2 - 2 * log( -theta[i] + 2 ) );

        break;

    default:
        // Code to execute if <variable> does not equal the value following any of the cases
        // Send a break message.
        break;
    }



}


template <bool propto__, bool jacobian__,
          typename Tu__, typename Tv__,typename Tt__, typename T__>
void bicop_log_latent(int i,
                   const std::vector<int>& copula_type,
                   Eigen::Matrix<Tu__,Eigen::Dynamic,1>& u,
                   Eigen::Matrix<Tv__,Eigen::Dynamic,1>&  v,
                   Eigen::Matrix<double,Eigen::Dynamic,1>& theta,
                   Eigen::Matrix<double,Eigen::Dynamic,1>& theta2,
                   T__& lp__,
                   stan::math::accumulator<T__>& lp_accum__,
                   stan::io::reader<T__>& in__
) {
    static const char* function("vifcopula::bicop_log_latent");
    using stan::math::square;
    using stan::math::log;

    switch ( copula_type[i] ) {
    case 0:
        // Independence copula
        lp_accum__.add(bicop_independence_log<propto__>(u,v));
        break;
    case 1:
        // Gaussian copula
        lp_accum__.add(bicop_normal_log<propto__>(u,v,theta[i]));
        break;
    case 2:
        // Student copula
        lp_accum__.add(bicop_student_log<propto__>(u,v,theta[i],theta2[i]));
        break;
    case 3:
        // Clayton copula
        lp_accum__.add(bicop_clayton_log<propto__>(u,v,theta[i]));
        break;
    case 4:
        // Gumbel copula
        lp_accum__.add(bicop_gumbel_log<propto__>(u,v,theta[i]));
        break;
    case 5:
        // Frank copula
        lp_accum__.add(bicop_frank_log<propto__>(u,v,theta[i]));
        break;
    case 6:
        // Joe copula
        lp_accum__.add(bicop_joe_log<propto__>(u,v,theta[i]));
        break;

    case 13:
        // survival Clayton copula
        lp_accum__.add(bicop_survival_clayton_log<propto__>(u,v,theta[i]));
        break;
    case 14:
        // survival Gumbel copula
        lp_accum__.add(bicop_survival_gumbel_log<propto__>(u,v,theta[i]));
        break;
    case 16:
        // survival Joe copula
        lp_accum__.add(bicop_survival_joe_log<propto__>(u,v,theta[i]));
        break;

    case 21:
        // Negative Guassian Copula
        lp_accum__.add(bicop_normal_log<propto__>(u,v,theta[i]));
        break;
    case 22:
        // Student copula
        lp_accum__.add(bicop_student_log<propto__>(u,v,theta[i],theta2[i]));
        break;

    case 23:
        // rotated 90 degree Clayton copula
        lp_accum__.add(bicop_r90_clayton_log<propto__>(u,v,theta[i]));
        break;
    case 24:
        // rotated 90 degree Gumbel copula
        lp_accum__.add(bicop_r90_gumbel_log<propto__>(u,v,theta[i]));
        break;

    case 25:
        // Frank copula
        lp_accum__.add(bicop_frank_log<propto__>(u,v,theta[i]));
        break;

    case 26:
        // rotated 90 degree Joe copula
        lp_accum__.add(bicop_r90_joe_log<propto__>(u,v,theta[i]));
        break;

    case 33:
        // rotated 270 degree Clayton copula
        lp_accum__.add(bicop_r270_clayton_log<propto__>(u,v,theta[i]));
        break;
    case 34:
        // rotated 270 degree Gumbel copula
        lp_accum__.add(bicop_r270_gumbel_log<propto__>(u,v,theta[i]));
        break;
    case 36:
        // rotated 270 degree Joe copula
        lp_accum__.add(bicop_r270_joe_log<propto__>(u,v,theta[i]));
        break;

    default:
        // Code to execute if <variable> does not equal the value following any of the cases
        // Send a break message.
        break;
    }



}


double bicop_log_double( const int copula_type,
                   double u,
                   double  v,
                   double theta,
                   double  theta2 = 0
) {
    static const char* function("vifcopula::bicop_log_double");
    double log_bicop = 0;
    switch ( copula_type ) {
    case 0:
        // Independence copula
        log_bicop = bicop_independence_log<FALSE>(u,v);
        break;
    case 1:
        // Gaussian copula
        log_bicop = bicop_normal_log<FALSE>(u,v,theta);
        break;
    case 2:
        // Student copula
        log_bicop = bicop_student_log<FALSE>(u,v,theta,theta2);
        break;
    case 3:
        // Clayon copula
        log_bicop = bicop_clayton_log<FALSE>(u,v,theta);
        break;
    case 4:
        // Gumbel copula
        log_bicop = bicop_gumbel_log<FALSE>(u,v,theta);
        break;
    case 5:
        // Frank copula
        log_bicop = bicop_frank_log<FALSE>(u,v,theta);
        break;
    case 6:
        // Joe copula
        log_bicop = bicop_joe_log<FALSE>(u,v,theta);
        break;

    case 13:
        // survival Clayon copula
        log_bicop = bicop_survival_clayton_log<FALSE>(u,v,theta);
        break;
    case 14:
        // survival Gumbel copula
        log_bicop = bicop_survival_gumbel_log<FALSE>(u,v,theta);
        break;
    case 16:
        // survival Joe copula
        log_bicop = bicop_survival_joe_log<FALSE>(u,v,theta);
        break;

    case 21:
        // Gaussian neg corr
        log_bicop = bicop_normal_log<FALSE>(u,v,theta);
        break;
    case 22:
        // Student copula neg corr

        log_bicop = bicop_student_log<FALSE>(u,v,theta,theta2);
        break;

    case 23:
        // rotated 90 degree Clayon copula
        log_bicop = bicop_r90_clayton_log<FALSE>(u,v,theta);
        break;
    case 24:
        // rotated 90 degree Gumbel copula
        log_bicop = bicop_r90_gumbel_log<FALSE>(u,v,theta);
        break;

    case 25:
        // Frank copula
        log_bicop = bicop_frank_log<FALSE>(u,v,theta);
        break;

    case 26:
        // rotated 90 degree Joe copula
        log_bicop = bicop_r90_joe_log<FALSE>(u,v,theta);
        break;

    case 33:
        // rotated 270 degree Clayon copula
        log_bicop = bicop_r270_clayton_log<FALSE>(u,v,theta);
        break;
    case 34:
        // rotated 270 degree Gumbel copula
        log_bicop = bicop_r270_gumbel_log<FALSE>(u,v,theta);
        break;
    case 36:
        // rotated 270 degree Joe copula
        log_bicop = bicop_r270_joe_log<FALSE>(u,v,theta);
        break;

    default:
        // Code to execute if <variable> does not equal the value following any of the cases
        // Send a break message.
        break;
    }

    return log_bicop;

} // function

}




#endif // VIFCOPULA_DISTRIBUTION_BICOP_LOG_HPP


