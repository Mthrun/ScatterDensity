// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// c_inPSphere2D
IntegerVector c_inPSphere2D(NumericMatrix data, IntegerVector xBinNr, IntegerVector yBinNr, unsigned int nrXBins, unsigned int nrYBins, unsigned int nrData, double paretoRadius);
RcppExport SEXP _ScatterDensity_c_inPSphere2D(SEXP dataSEXP, SEXP xBinNrSEXP, SEXP yBinNrSEXP, SEXP nrXBinsSEXP, SEXP nrYBinsSEXP, SEXP nrDataSEXP, SEXP paretoRadiusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type xBinNr(xBinNrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type yBinNr(yBinNrSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nrXBins(nrXBinsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nrYBins(nrYBinsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nrData(nrDataSEXP);
    Rcpp::traits::input_parameter< double >::type paretoRadius(paretoRadiusSEXP);
    rcpp_result_gen = Rcpp::wrap(c_inPSphere2D(data, xBinNr, yBinNr, nrXBins, nrYBins, nrData, paretoRadius));
    return rcpp_result_gen;
END_RCPP
}
// c_quantile
Rcpp::NumericVector c_quantile(Rcpp::NumericVector x, Rcpp::NumericVector probs);
RcppExport SEXP _ScatterDensity_c_quantile(SEXP xSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(c_quantile(x, probs));
    return rcpp_result_gen;
END_RCPP
}
// insidecpp
bool insidecpp(const NumericVector& xy, const int n1, const int n2, const NumericMatrix& poly);
RcppExport SEXP _ScatterDensity_insidecpp(SEXP xySEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP polySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type xy(xySEXP);
    Rcpp::traits::input_parameter< const int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< const int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type poly(polySEXP);
    rcpp_result_gen = Rcpp::wrap(insidecpp(xy, n1, n2, poly));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ScatterDensity_c_inPSphere2D", (DL_FUNC) &_ScatterDensity_c_inPSphere2D, 7},
    {"_ScatterDensity_c_quantile", (DL_FUNC) &_ScatterDensity_c_quantile, 2},
    {"_ScatterDensity_insidecpp", (DL_FUNC) &_ScatterDensity_insidecpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_ScatterDensity(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
