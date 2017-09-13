#ifndef VIFCOPULA_EXTRA_VINECOPULA_HPP
#define VIFCOPULA_EXTRA_VINECOPULA_HPP

#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::interfaces(r,cpp)]]

void (*Hfunc2) (int* family,int* n,double* v,double* u,double* theta,double* nu,double* out);

void (*diffhfunc_rho_tCopula) (double* u, double* v, int* n, double* param, int* copula, double* out);      // par
void (*diffhfunc_mod) (double* u, double* v, int* n, double* param, int* copula, double* out);              // par
void (*diffhfunc_nu_tCopula_new) (double* u, double* v, int* n, double* param, int* copula, double* out);   // par2
void (*diffhfunc_v_mod) (double* u, double* v, int* n, double* param, int* copula, double* out);            // u2
void (*difflPDF_nu_tCopula_new) (double* u, double* v, int* n, double* param, int* copula, double* out);    // deriv par2


extern "C" void R_init_vifcopula_extra(DllInfo *dll) {
    Hfunc2 = (void (*) (int* ,int* ,double* ,double* ,double* ,double* ,double* )) R_GetCCallable("VineCopula", "Hfunc2");

    diffhfunc_rho_tCopula = (void (*) (double* , double* , int* , double* , int* , double* )) R_GetCCallable("VineCopula", "diffhfunc_rho_tCopula");
    diffhfunc_mod = (void (*) (double* , double* , int* , double* , int* , double* )) R_GetCCallable("VineCopula", "diffhfunc_mod");

    diffhfunc_nu_tCopula_new = (void (*) (double* , double* , int* , double* , int* , double* )) R_GetCCallable("VineCopula", "diffhfunc_nu_tCopula_new");

    diffhfunc_v_mod = (void (*) (double* , double* , int* , double* , int* , double* )) R_GetCCallable("VineCopula", "diffhfunc_v_mod");

    difflPDF_nu_tCopula_new = (void (*) (double* , double* , int* , double* , int* , double* )) R_GetCCallable("VineCopula", "difflPDF_nu_tCopula_new");
}

#endif // VIFCOPULA_EXTRA_VINECOPULA_HPP
