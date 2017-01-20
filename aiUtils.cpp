#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

#define NA 9

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rowMedian_cpp(arma::mat &x){
  int nrow = x.n_rows ;
  int ncol = x.n_cols ;
  NumericVector out(nrow);
  for (int i = 0; i < nrow; i++){ 
    int len_row_i = 0;
    for (int j = 0; j < ncol; j++){
      if (!NumericVector::is_na(x.at(i, j)) && (x.at(i, j) != NA)){
        len_row_i ++;  
      }      
    }
    NumericVector row_i(len_row_i);
    int count = 0;
    for (int j = 0; j < ncol; j++){
      if (!NumericVector::is_na(x.at(i, j)) && (x.at(i, j) != NA)){
        row_i[count] = x.at(i, j);
        count ++;  
      }      
    }
    out[i] = Rcpp::median(row_i);
  }
  return out;
}

// [[Rcpp::export]]
arma::mat impute_geno(arma::mat &x){
  int nrow = x.n_rows ;
  int ncol = x.n_cols ;
  arma::mat y = x;
  NumericVector rowmed = rowMedian_cpp(x);
  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      if (NumericVector::is_na(x.at(i, j)) || (x.at(i, j) == NA)){
        y.at(i, j) = rowmed[i];  
      }
    }
  }
  return y;
}