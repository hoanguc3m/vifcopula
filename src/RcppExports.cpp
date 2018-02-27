// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/vifcopula.h"
#include <RcppEigen.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// vifcop
List vifcop(SEXP data_, SEXP init_, SEXP other_);
static SEXP vifcopula_vifcop_try(SEXP data_SEXP, SEXP init_SEXP, SEXP other_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type data_(data_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type init_(init_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type other_(other_SEXP);
    rcpp_result_gen = Rcpp::wrap(vifcop(data_, init_, other_));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP vifcopula_vifcop(SEXP data_SEXP, SEXP init_SEXP, SEXP other_SEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(vifcopula_vifcop_try(data_SEXP, init_SEXP, other_SEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// hmcfcop
List hmcfcop(SEXP data_, SEXP init_, SEXP other_);
static SEXP vifcopula_hmcfcop_try(SEXP data_SEXP, SEXP init_SEXP, SEXP other_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type data_(data_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type init_(init_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type other_(other_SEXP);
    rcpp_result_gen = Rcpp::wrap(hmcfcop(data_, init_, other_));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP vifcopula_hmcfcop(SEXP data_SEXP, SEXP init_SEXP, SEXP other_SEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(vifcopula_hmcfcop_try(data_SEXP, init_SEXP, other_SEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int vifcopula_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("List(*vifcop)(SEXP,SEXP,SEXP)");
        signatures.insert("List(*hmcfcop)(SEXP,SEXP,SEXP)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP vifcopula_RcppExport_registerCCallable() { 
    R_RegisterCCallable("vifcopula", "vifcopula_vifcop", (DL_FUNC)vifcopula_vifcop_try);
    R_RegisterCCallable("vifcopula", "vifcopula_hmcfcop", (DL_FUNC)vifcopula_hmcfcop_try);
    R_RegisterCCallable("vifcopula", "vifcopula_RcppExport_validate", (DL_FUNC)vifcopula_RcppExport_validate);
    return R_NilValue;
}
