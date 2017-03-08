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
    ber = R::rbinom(2, p);
    if (ber == 1){
      g = 0;
      jumps[i] = 1;
    }
  }
  return(jumps);
}

//' @export
//' 
// [[Rcpp::export]]
LogicalVector ancestry_chunks(const NumericVector &jumps){
  int n = jumps.size();  
  LogicalVector chunks(n);
  bool fbool = true;
  for (int i = 0; i < n; i++){
    if (jumps[i] == 1){
      fbool = !fbool;
    }
    chunks[i] = fbool;
  }
  return(chunks);
}

//' @export
//' 
// [[Rcpp::export]]
NumericVector generate_hybrid_cpp(const NumericMatrix &H1, const NumericMatrix &H2, const double alpha, const LogicalVector chunks){
  int n = chunks.size();
  NumericVector haplotype_1(n);
  NumericVector haplotype_2(n);
  double p = alpha;
  double nbino = R::rbinom(1.0, p);
  for (int i = 0; i < n; i++){
    if (chunks[i]){
      p = 1 - p;
      nbino = R::rbinom(1.0, p);
    }
  }
  return(haplotype_1);
}

//' @export
//' 
// [[Rcpp::export]]
NumericMatrix haplo_to_geno(const NumericMatrix &H){
  int nSNP = H.nrow();
  int nIND = H.ncol() / 2;
  NumericMatrix G(nSNP, nIND);
  for (int j = 0; j < nIND; j++){
    for (int i = 0; i < nSNP; i++){
      G(i, j) = H(i, 2 * j) + H(i, 2 * j + 1); 
    }
  }
  return G;
}
