// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cpt_uni_modf
int cpt_uni_modf(NumericVector chain);
RcppExport SEXP _cptmcmc_cpt_uni_modf(SEXP chainSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type chain(chainSEXP);
    rcpp_result_gen = Rcpp::wrap(cpt_uni_modf(chain));
    return rcpp_result_gen;
END_RCPP
}
// cpt_multi_modf
NumericVector cpt_multi_modf(NumericMatrix chain);
RcppExport SEXP _cptmcmc_cpt_multi_modf(SEXP chainSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type chain(chainSEXP);
    rcpp_result_gen = Rcpp::wrap(cpt_multi_modf(chain));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cptmcmc_cpt_uni_modf", (DL_FUNC) &_cptmcmc_cpt_uni_modf, 1},
    {"_cptmcmc_cpt_multi_modf", (DL_FUNC) &_cptmcmc_cpt_multi_modf, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_cptmcmc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
