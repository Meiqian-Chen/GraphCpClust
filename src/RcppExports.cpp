// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// GridSearch1
Rcpp::List GridSearch1(const NumericVector x, const NumericVector y, const float beta00, const float beta01, const float beta10, const float beta11, const int n);
RcppExport SEXP _GraphCpClust_GridSearch1(SEXP xSEXP, SEXP ySEXP, SEXP beta00SEXP, SEXP beta01SEXP, SEXP beta10SEXP, SEXP beta11SEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const float >::type beta00(beta00SEXP);
    Rcpp::traits::input_parameter< const float >::type beta01(beta01SEXP);
    Rcpp::traits::input_parameter< const float >::type beta10(beta10SEXP);
    Rcpp::traits::input_parameter< const float >::type beta11(beta11SEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(GridSearch1(x, y, beta00, beta01, beta10, beta11, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GraphCpClust_GridSearch1", (DL_FUNC) &_GraphCpClust_GridSearch1, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_GraphCpClust(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}