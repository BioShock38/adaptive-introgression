#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

#define NA 9

using namespace Rcpp;


// [[Rcpp::export]]
arma::mat cmpt_global_pca(arma::mat geno, arma::mat V, arma::vec sigma){
  int nIND = geno.n_rows;
  int nSNP = V.n_rows;
  int K = V.n_cols;
  arma::mat u(nIND, K);
  u.zeros();
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      for (int i = 0; i < nSNP; i++){
        u(j, k) += geno(j, i) * V(i, k) / sigma[k];
      }
    }
  }
  return(u);
}

// [[Rcpp::export]]
arma::mat cmpt_local_pca(arma::mat geno, arma::mat V, arma::vec sigma, int beg, int end){
  int nIND = geno.n_rows;
  int nSNP = V.n_rows;
  int K = V.n_cols;
  arma::mat uloc(nIND, K);
  uloc.zeros();
  double cst = (double) nSNP;
  cst /=  (double) end - beg;
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      for (int i = beg; i < end; i++){
        uloc(j, k) += geno(j, i) * V(i, k) * cst / sigma[k];
      }
    }
  }
  return(uloc);
}

// [[Rcpp::export]]
int get_nb_ind(arma::vec lab, int anc){
  int n = lab.n_elem;
  int c = 0;
  for (int i = 0; i < n; i++){
    if (lab[i] == anc){
      c += 1;  
    }
  }
  return(c);
}

// [[Rcpp::export]]
Rcpp::List cmpt_centroids(arma::mat u, arma::vec lab, int anc1, int anc2){
  int nIND = u.n_rows;
  int K = u.n_cols;
  arma::vec m1(K);
  m1.zeros();
  arma::vec m2(K);
  m2.zeros();
  int c1 = get_nb_ind(lab, anc1);
  int c2 = get_nb_ind(lab, anc2);
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      if (lab[j] == anc1){
        m1[k] += u(j, k) / c1;
      } else if (lab[j] == anc2){
        m2[k] += u(j, k) / c2;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("m1") = m1,
                            Rcpp::Named("m2") = m2);
}

// [[Rcpp::export]]
arma::mat rescale_local_pca(arma::mat &u, arma::vec &s, arma::vec &dep_glob, arma::vec &dep_loc){
  int nIND = u.n_rows;
  int K = u.n_cols;
  arma::mat usc(nIND, K);
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      usc(j, k) = u(j, k) * s[k];
      usc(j, k) = usc(j, k) + dep_glob[k] - dep_loc[k];
    }
  }
  return(usc);
}



// [[Rcpp::export]]
double cmpt_window_stat(arma::mat &geno, 
                           arma::mat &V, 
                           arma::vec &sigma, 
                           arma::mat &uglob, 
                           int beg, 
                           int end, 
                           int direction, 
                           arma::vec lab, 
                           int adm, 
                           int axis){
  arma::mat uloc = cmpt_local_pca(geno, V, sigma, beg, end);
  int nADM = uglob.n_rows; 
  double stat = 0;
  if (direction == 1){
    for (int j = 0; j < nADM; j++){
      if ((lab[j] == adm) && (uloc(j, axis) - uglob(j, axis)) > 0){
        stat += (uloc(j, axis) - uglob(j, axis)) * (uloc(j, axis) - uglob(j, axis));
      }
    }
  } else if (direction == (-1)){
    for (int j = 0; j < nADM; j++){
      if ((lab[j] == adm) && (uloc(j, axis) - uglob(j, axis)) < 0){
        stat += (uloc(j, axis) - uglob(j, axis)) * (uloc(j, axis) - uglob(j, axis));
      }
    }
  }
  return(stat);
}
// [[Rcpp::export]]
arma::vec cmpt_all_stat(arma::mat &geno, 
                        arma::mat &V, 
                        arma::vec &sigma, 

                        int direction, 
                        arma::vec lab, 
                        int adm, 
                        int axis){
  int nSNP = geno.n_cols;
  for (int i = 0; i < (nSNP - window_size); i++){
    
  }

}

