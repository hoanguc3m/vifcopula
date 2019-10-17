#ifndef VIFCOPULA_HFUNC_STAN_HPP
#define VIFCOPULA_HFUNC_STAN_HPP

#include <stan/math/rev/core.hpp>

// [[Rcpp::interfaces(r,cpp)]]


void (*Hfunc2) (int* family,int* n,double* v,double* u,double* theta,double* nu,double* out);

void (*diffhfunc_rho_tCopula) (double* u, double* v, int* n, double* param, int* copula, double* out);      // par
void (*diffhfunc_mod) (double* u, double* v, int* n, double* param, int* copula, double* out);              // par
void (*diffhfunc_nu_tCopula_new) (double* u, double* v, int* n, double* param, int* copula, double* out);   // par2
void (*diffhfunc_v_mod) (double* u, double* v, int* n, double* param, int* copula, double* out);            // u2


namespace vifcopula
{
using namespace stan::math;
using namespace vifcopula;

template <typename T> int sgn(T val) {
    return (val < T(20)) - (T(20) < val);
}

double diffhfuncBB1_v(int family,
                      double u,
                      double v,
                      double theta,
                      double delta,
                      double hfunc2_val){
    double hfunc2_dv = 0;
    double sign = 1.;

    if (family == 17){
        u = 1. - u;
        v = 1. - v;
    }

    if (family == 27){
        u = 1. - u;
        theta = - theta;
        delta = - delta;
        sign = -1.;
    }
    if (family == 37){
        v = 1. - v;
        theta = - theta;
        delta = - delta;
        sign = -1.;
    }

    // if (family == 7){
    double x = pow(u, -theta) - 1;
    double y = pow(v, -theta) - 1;              // A3
    double A3 = y;
    double A2 = pow(x, delta)+pow(y, delta);    // A2
    double A1 = 1 + pow(A2,1./delta);           // A1

    // hfunc = pow(A1, -1/theta - 1) * pow(A2, 1/delta - 1) * pow(A3, delta-1) * pow(u2, -theta-1);
    double B4 = pow(v, -theta-1);
    double dydv = -theta * B4;
    double B3 = pow(y, delta-1);
    double B2 = pow(A2, 1./delta - 1);
    double B1 = pow(A1, -1./theta - 1);




    hfunc2_dv = -(theta+1) * (B4 / v) * B1 * B2 * B3 +         \
        (delta-1) * (B3/A3) * dydv * B1 * B2 * B4 +     \
        (1- delta) * (B2 / A2) * B3 * dydv * B1 * B3 * B4 + \
        - (1./theta +1) * (B1 / A1) * B2 * B3 * dydv * B2 * B3 * B4;
        // }

        return(sign*hfunc2_dv);

}

double diffhfuncBB1_theta(int family,
                          double u,
                          double v,
                          double theta,
                          double theta2,
                          double hfunc2_val){
    int rotated = sgn(family);

    double eps = 0.001 * rotated;
    double hfunc2_dtheta2 = 0;

    theta = theta + eps;
    double hfunc2_val_eps;
    int nn = 1;

    Hfunc2(&family, &nn, &u, &v, &theta, &theta2, &hfunc2_val_eps);
    hfunc2_dtheta2 = rotated * (hfunc2_val_eps - hfunc2_val) / eps;
    return(hfunc2_dtheta2);
}

double diffhfuncBB1_delta(int family,
                          double u,
                          double v,
                          double theta,
                          double theta2,
                          double hfunc2_val){
    int rotated = sgn(family);

    double eps = 0.001 * rotated;
    double hfunc2_dtheta2;
    int nn = 1;
    theta2 = theta2 + eps;
    double hfunc2_val_eps;
    Hfunc2(&family, &nn, &u, &v, &theta, &theta2, &hfunc2_val_eps);
    hfunc2_dtheta2 = rotated * (hfunc2_val_eps - hfunc2_val) / eps;
    return(hfunc2_dtheta2);
}

// hfunc_trans with 3 vars arguement
stan::math::var hfunc_trans (int family,
                             double u,
                             const stan::math::var& v,
                             const stan::math::var& theta,
                             const stan::math::var& theta2){
    double v_val = v.val();
    double theta_val = theta.val();
    double theta2_val = theta2.val();
    double hfunc2_val;
    double hfunc2_dv = 0.;
    double hfunc2_dtheta = 0.;
    double hfunc2_dtheta2 = 0.;



    int nn = 1;

    Hfunc2(&family, &nn, &u, &v_val, &theta_val, &theta2_val, &hfunc2_val);

    double *param = new double[2];
    param[0] = theta_val;
    param[1] = theta2_val;

    // Derivative of Hfunc2 of wrt v
    if (family % 10 == 7){
        hfunc2_dv = diffhfuncBB1_v(family, u, v_val, theta_val, theta2_val, hfunc2_val);
    } else {
        diffhfunc_v_mod(&u, &v_val, &nn, param, &family, &hfunc2_dv);
    }

    // Derivative of Hfunc2 of Student copula wrt rho
    if (family % 10 == 2){
        diffhfunc_rho_tCopula(&u, &v_val, &nn, param, &family, &hfunc2_dtheta);
    } else if (family % 10 == 7){
        // Derivative of Hfunc2 of BB1 copula wrt theta
        hfunc2_dtheta = diffhfuncBB1_theta(family, u, v_val, theta_val, theta2_val, hfunc2_val);
    } else {
        diffhfunc_mod(&u, &v_val, &nn, param, &family, &hfunc2_dtheta);
    }

    // Derivative wrt to the par2
    // Derivative of Hfunc2 of Student copula wrt nu
    if (family % 10 == 2){
        diffhfunc_nu_tCopula_new(&u, &v_val, &nn, param, &family, &hfunc2_dtheta2);
    }
    // Derivative of Hfunc2 of BB1 copula wrt delta
    if (family % 10 == 7){
        hfunc2_dtheta2 = diffhfuncBB1_delta(family, u, v_val, theta_val, theta2_val, hfunc2_val);
    }


    delete param;
    return (new precomp_vvv_vari(hfunc2_val, v.vi_, theta.vi_, theta2.vi_, hfunc2_dv, hfunc2_dtheta, hfunc2_dtheta2) );

}

// hfunc_trans with 2 vars arguement
stan::math::var hfunc_trans (int family,
                             double u,
                             const stan::math::var& v,
                             const stan::math::var& theta,
                             double theta2){
    double v_val = v.val();
    double theta_val = theta.val();
    double hfunc2_val;
    double hfunc2_dv = 0.;
    double hfunc2_dtheta = 0.;



    int nn = 1;

    Hfunc2(&family, &nn, &u, &v_val, &theta_val, &theta2, &hfunc2_val);

    double *param = new double[2];
    param[0] = theta_val;
    param[1] = theta2;

    // Derivative of Hfunc2 of wrt v
    if (family % 10 == 7){
        hfunc2_dv = diffhfuncBB1_v(family, u, v_val, theta_val, theta2, hfunc2_val);
    } else {
        diffhfunc_v_mod(&u, &v_val, &nn, param, &family, &hfunc2_dv);
    }

    // Derivative of Hfunc2 of Student copula wrt rho
    if (family % 10 == 2){
        diffhfunc_rho_tCopula(&u, &v_val, &nn, param, &family, &hfunc2_dtheta);
    } else if (family % 10 == 7){
        // Derivative of Hfunc2 of BB1 copula wrt theta
        hfunc2_dtheta = diffhfuncBB1_theta(family, u, v_val, theta_val, theta2, hfunc2_val);
    } else {
        diffhfunc_mod(&u, &v_val, &nn, param, &family, &hfunc2_dtheta);
    }

    delete param;
    return (new precomp_vv_vari(hfunc2_val, v.vi_, theta.vi_, hfunc2_dv, hfunc2_dtheta) );
}
// hfunc_trans with 2 vars arguement
stan::math::var hfunc_trans (int family,
                             double u,
                             double v,
                             const stan::math::var& theta,
                             const stan::math::var& theta2){
    double theta_val = theta.val();
    double theta2_val = theta2.val();
    double hfunc2_val;
    double hfunc2_dtheta = 0.;
    double hfunc2_dtheta2 = 0.;


    int nn = 1;

    Hfunc2(&family, &nn, &u, &v, &theta_val, &theta2_val, &hfunc2_val);

    double *param = new double[2];
    param[0] = theta_val;
    param[1] = theta2_val;

    // Derivative of Hfunc2 of Student copula wrt rho
    if (family % 10 == 2){
        diffhfunc_rho_tCopula(&u, &v, &nn, param, &family, &hfunc2_dtheta);
    } else if (family % 10 == 7){
        // Derivative of Hfunc2 of BB1 copula wrt theta
        hfunc2_dtheta = diffhfuncBB1_theta(family, u, v, theta_val, theta2_val, hfunc2_val);
    } else {
        diffhfunc_mod(&u, &v, &nn, param, &family, &hfunc2_dtheta);
    }

    // Derivative wrt to the par2
    // Derivative of Hfunc2 of Student copula wrt nu
    if (family % 10 == 2){
        diffhfunc_nu_tCopula_new(&u, &v, &nn, param, &family, &hfunc2_dtheta2);
    }
    // Derivative of Hfunc2 of BB1 copula wrt delta
    if (family % 10 == 7){
        hfunc2_dtheta2 = diffhfuncBB1_delta(family, u, v, theta_val, theta2_val, hfunc2_val);
    }



    delete param;
    return (new precomp_vv_vari(hfunc2_val, theta.vi_, theta2.vi_, hfunc2_dtheta, hfunc2_dtheta2) );
}



// hfunc_trans with 1 vars arguement - v vari
stan::math::var hfunc_trans (int family,
                             double u,
                             const stan::math::var& v,
                             double theta,
                             double theta2){
    double v_val = v.val();
    double hfunc2_val;
    double hfunc2_dv = 0.;


    int nn = 1;

    Hfunc2(&family, &nn, &u, &v_val, &theta, &theta2, &hfunc2_val);

    double *param = new double[2];
    param[0] = theta;
    param[1] = theta2;

    // Derivative of Hfunc2 of wrt v
    if (family % 10 == 7){
        hfunc2_dv = diffhfuncBB1_v(family, u, v_val, theta, theta2, hfunc2_val);
    } else {
        diffhfunc_v_mod(&u, &v_val, &nn, param, &family, &hfunc2_dv);
    }

    delete param;
    return (new precomp_v_vari(hfunc2_val, v.vi_, hfunc2_dv) );

}


// hfunc_trans with 1 vars arguement - theta vari
stan::math::var hfunc_trans (int family,
                             double u,
                             double v,
                             const stan::math::var& theta,
                             double theta2){
    double theta_val = theta.val();
    double hfunc2_val;
    double hfunc2_dtheta = 0.;


    int nn = 1;

    Hfunc2(&family, &nn, &u, &v, &theta_val, &theta2, &hfunc2_val);

    double *param = new double[2];
    param[0] = theta_val;
    param[1] = theta2;

    // Derivative of Hfunc2 of Student copula wrt rho
    if (family % 10 == 2){
        diffhfunc_rho_tCopula(&u, &v, &nn, param, &family, &hfunc2_dtheta);
    } else if (family % 10 == 7){
        // Derivative of Hfunc2 of BB1 copula wrt theta
        hfunc2_dtheta = diffhfuncBB1_theta(family, u, v, theta_val, theta2, hfunc2_val);
    } else {
        diffhfunc_mod(&u, &v, &nn, param, &family, &hfunc2_dtheta);
    }


    delete param;
    return (new precomp_v_vari(hfunc2_val, theta.vi_, hfunc2_dtheta) );
}


// hfunc_trans with 1 vars arguement - theta2 vari
stan::math::var hfunc_trans (int family,
                             double u,
                             double v,
                             double theta,
                             const stan::math::var& theta2){
    double theta2_val = theta2.val();
    double hfunc2_val;
    double hfunc2_dtheta2 = 0.;


    int nn = 1;

    Hfunc2(&family, &nn, &u, &v, &theta, &theta2_val, &hfunc2_val);

    double *param = new double[2];
    param[0] = theta;
    param[1] = theta2_val;

    // Derivative of Hfunc2 of Student copula wrt nu
    if (family % 10 == 2){
        diffhfunc_nu_tCopula_new(&u, &v, &nn, param, &family, &hfunc2_dtheta2);
    }
    // Derivative of Hfunc2 of BB1 copula wrt delta
    if (family % 10 == 7){
        hfunc2_dtheta2 = diffhfuncBB1_delta(family, u, v, theta, theta2_val, hfunc2_val);
    }

    delete param;
    return (new precomp_v_vari(hfunc2_val, theta2.vi_, hfunc2_dtheta2) );

}

// hfunc_trans with 0 vars arguement
double hfunc_trans (int family,
                    double u,
                    double v,
                    double theta,
                    double theta2){
    double hfunc2_val;
    int nn = 1;

    Hfunc2(&family, &nn, &u, &v, &theta, &theta2, &hfunc2_val);

    return hfunc2_val;

}

}       // end name_sapce
#endif // VIFCOPULA_HFUNC_HPP
