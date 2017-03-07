#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;


Rcpp::List multipleLm(arma::mat &Y, arma::mat &X) {
  
  int n = X.n_rows , k = X.n_cols;
  int p = Y.n_cols;
  
  arma::colvec coef = arma::solve(X, y);      // fit model y ~ X
  arma::colvec resid = y - X*coef;            // residuals
  
  double sig2 = arma::as_scalar( arma::trans(resid)*resid/(n-k) );
  // std.error of estimate
  arma::colvec stderrest = arma::sqrt( sig2 * arma::diagvec( arma::inv(arma::trans(X)*X)) );
  
  return Rcpp::List::create(
    Rcpp::Named("coefficients") = coef,
    Rcpp::Named("stderr")       = stderrest
  ) ;
  
}