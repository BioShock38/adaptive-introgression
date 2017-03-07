#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

//' @export
//' 
// [[Rcpp::export]]
NumericVector jumps_from_map(const NumericVector &map, const double &lambda){
  int n = map.size();  
  NumericVector jumps(n);
  double g = 0;
  double p = 0;
  double ber = 0;
  for (int i = 0; i < n; i++){
    g += map[i];
    p = R::pexp(g, lambda, 1, 0);
    ber = R::rbinom(1, p);
    if (ber == 1){
      g = 0;
      jumps[i] = 1;
    }
  }
  return(jumps);
}
