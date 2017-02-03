#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

#define NA 9

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat compute_stat(arma::mat &geno, arma::vec &scores, arma::vec &loadings, double sigma, int window_size, int direction){
  int nSNP = geno.n_cols;
  int nIND = geno.n_rows;
  arma::vec stat(nSNP);
  stat.zeros();
  arma::vec tmp(nIND);
  tmp.zeros();
  
  double sum_loadings = 0;
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < window_size; k++){
      tmp[j] += geno.at(j, k) * loadings[k] / sigma;
      sum_loadings += loadings[k] * loadings[k];
    }
    if (direction == 1){
      if ((tmp[j] - scores[j]) > 0){
        stat[0] += ((tmp[j] / sum_loadings) - scores[j]) * ((tmp[j] / sum_loadings) - scores[j]); 
      }
    } else if (direction == -1){
      if ((tmp[j] - scores[j]) < 0){
        stat[0] += ((tmp[j] / sum_loadings) - scores[j]) * ((tmp[j] / sum_loadings) - scores[j]); 
      }
    } else if (direction == 0){
      stat[0] += ((tmp[j] / sum_loadings) - scores[j]); 
    }
  }
  
  for (int i = 1; i < (nSNP - window_size); i++){
    for (int j = 0; j < nIND; j++){
      int i_beg = i - 1;
      int i_end = i + window_size;
      tmp[j] -= geno.at(j, i_beg) * loadings[i_beg] / sigma;
      sum_loadings -= loadings[i_beg] * loadings[i_beg]; 
      tmp[j] += geno.at(j, i_end) * loadings[i_end] / sigma;
      sum_loadings += loadings[i_end] * loadings[i_end];
      if (direction == 1){
        if ((tmp[j] - scores[j]) > 0){
          stat[i] += ((tmp[j] / sum_loadings) - scores[j]) * ((tmp[j] / sum_loadings) - scores[j]); 
        }
      } else if (direction == -1){
        if ((tmp[j] - scores[j]) < 0){
          stat[i] += ((tmp[j] / sum_loadings) - scores[j]) * ((tmp[j] / sum_loadings) - scores[j]); 
        }
      } else if (direction == 0){
        stat[i] += ((tmp[j] / sum_loadings) - scores[j]); 
      }
    }
    if (direction == 0){
      stat[i] *= stat[i];
    }
  }
  for (int i = (nSNP - window_size); i < nSNP; i++){
    stat[i] = stat[nSNP - window_size - 1];
  }
  return(stat);
}


// [[Rcpp::export]]
arma::mat compute_stat_kernel(arma::mat &geno, arma::vec &scores, arma::vec &loadings, double sigma, int direction, arma::mat &kernel){
  int nSNP = geno.n_cols;
  int nIND = geno.n_rows;
  arma::vec stat(nSNP);
  stat.zeros();
  int window_size = 25000;
  arma::vec tmp(nIND);
  arma::vec mask(nIND);
  arma::vec lambda;
  tmp.zeros();
  mask.zeros();
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < window_size; k++){
      tmp[j] += geno.at(j, k) * loadings[k] / sigma;
    }
    tmp[j] *= double (nSNP / window_size);
    if (direction == 1){
      if ((tmp[j] - scores[j]) > 0){
        mask[j] = tmp[j];
      } else {
        mask[j] = 0.0;
      }
    } else if (direction == 0){
      if ((tmp[j] - scores[j]) < 0){
        mask[j] = tmp[j];
      } else {
        mask[j] = 0.0;
      }
    }
  }
  lambda = arma::solve(kernel, mask);
  stat[0] = dot(lambda, kernel * lambda);
  for (int i = 1; i < (nSNP - window_size); i++){
    mask.zeros();
    for (int j = 0; j < nIND; j++){
      int i_beg = i - 1;
      int i_end = i + window_size;
      tmp[j] -= geno.at(j, i_beg) * loadings[i_beg] * nSNP/ (window_size * sigma);
      tmp[j] += geno.at(j, i_end) * loadings[i_end] * nSNP/ (window_size * sigma);
      if (direction == 1){
        if ((tmp[j] - scores[j]) > 0){
          mask[j] = tmp[j];
        } else {
          mask[j] = 0.0;
        }
      } else if (direction == 0){
        if ((tmp[j] - scores[j]) < 0){
          mask[j] = tmp[j];
        } else {
          mask[j] = 0.0;
        }
      }
    }
    lambda = arma::solve(kernel, mask);
    stat[i] = dot(lambda, kernel * lambda);
  }
  return(stat);
}


// [[Rcpp::export]]
arma::vec compute_stat_procrusted(arma::mat &geno_adm, 
                                  arma::mat &geno_anc, 
                                  arma::vec &scores_adm, 
                                  arma::vec &scores_anc, 
                                  arma::vec &loadings, 
                                  double sigma, 
                                  int window_size, 
                                  int direction){
  int nSNP = geno_adm.n_cols;
  int nind_adm = geno_adm.n_rows;
  int nind_anc = geno_anc.n_rows;
  arma::vec stat(nSNP);
  stat.zeros();
  
  /* Projecting individuals using a smaller subset of markers may affect the scores */
  /* The loadings have to be recalculated */ 
  arma::vec local_loadings(window_size);
  arma::vec tmp(nind_adm);
  
  for (int i = 0; i < (nSNP - window_size); i ++){
    int i_beg = i;
    int i_end = i + window_size - 1;
    local_loadings.zeros();
    tmp.zeros();
    local_loadings = arma::solve(geno_anc.submat(0, i_beg, nind_anc - 1, i_end), scores_anc);
    for (int j = 0; j < nind_adm; j++){
      for (int k = 0; k < window_size; k++){
        tmp[j] += geno_adm.at(j, k) * local_loadings[k];
      }
      if (direction == 1){
        if ((tmp[j] - scores_adm[j]) > 0){
          stat[i] += (tmp[j] - scores_adm[j]) * (tmp[j] - scores_adm[j]);
        }
      } else if (direction == 0){
        if ((tmp[j] - scores_adm[j]) < 0){
          stat[i] += (tmp[j] - scores_adm[j]) * (tmp[j] - scores_adm[j]); 
        }
      }
    }
  }
  return(stat);
}


// [[Rcpp::export]]
arma::mat compute_stat_local_pca(arma::mat &geno, arma::vec &scores, arma::vec &loadings, double sigma, int window_size, int direction){
  int nSNP = geno.n_cols;
  int nIND = geno.n_rows;
  arma::vec stat(nSNP);
  stat.zeros();
  arma::vec tmp(nIND);
  tmp.zeros();
  int n_count = 0;
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < window_size; k++){
      tmp[j] += geno.at(j, k) * loadings[k] / sigma;
    }
    tmp[j] *= double (nSNP / window_size);
    if (direction == 1){
      if ((tmp[j] - scores[j]) > 0){
        stat[0] += (tmp[j] - scores[j]) * (tmp[j] - scores[j]); 
        n_count++;
      }
    } else if (direction == 0){
      if ((tmp[j] - scores[j]) < 0){
        stat[0] += (tmp[j] - scores[j]) * (tmp[j] - scores[j]);  
        n_count++;
      }
    }
  }
  //stat[0] /= sqrt(nIND + 1 - n_count);
  for (int i = 1; i < (nSNP - window_size); i++){
    int n_count = 0;
    for (int j = 0; j < nIND; j++){
      int i_beg = i - 1;
      int i_end = i + window_size;
      tmp[j] -= geno.at(j, i_beg) * loadings[i_beg] * nSNP/ (window_size * sigma);
      tmp[j] += geno.at(j, i_end) * loadings[i_end] * nSNP/ (window_size * sigma);
      if (direction == 1){
        if ((tmp[j] - scores[j]) > 0){
          stat[i] += (tmp[j] - scores[j]) * (tmp[j] - scores[j]);   
          n_count++;
        }
      } else if (direction == -1){
        if ((tmp[j] - scores[j]) < 0){
          stat[i] += (tmp[j] - scores[j]) * (tmp[j] - scores[j]);   
          n_count++;
        }
      } else if (direction == 0){
        stat[i] += (tmp[j] - scores[j]);   
        n_count++;
      }
    }
    if (direction == 0){
      stat[i] *= stat[i];
    }
    //stat[i] /= sqrt(nIND + 1 - n_count);
  }
  for (int i = (nSNP - window_size); i < nSNP; i++){
    stat[i] = stat[nSNP - window_size - 1];
  }
  return(stat);
}

// [[Rcpp::export]]
arma::mat compute_stat(arma::mat &geno, arma::vec &scores, arma::vec &loadings, double sigma, int window_size, int direction){
  int nSNP = geno.n_cols;
  int nIND = geno.n_rows;
  arma::vec stat(nSNP);
  stat.zeros();
  arma::vec tmp(nIND);
  tmp.zeros();
  int n_count = 0;
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < window_size; k++){
      tmp[j] += geno.at(j, k) * loadings[k] / sigma;
    }
    tmp[j] *= double (nSNP / window_size);
    if (direction == 1){
      if ((tmp[j] - scores[j]) > 0){
        stat[0] += (tmp[j] - scores[j]) * (tmp[j] - scores[j]); 
        n_count++;
      }
    } else if (direction == 0){
      if ((tmp[j] - scores[j]) < 0){
        stat[0] += (tmp[j] - scores[j]) * (tmp[j] - scores[j]);  
        n_count++;
      }
    }
  }
  //stat[0] /= sqrt(nIND + 1 - n_count);
  for (int i = 1; i < (nSNP - window_size); i++){
    int n_count = 0;
    for (int j = 0; j < nIND; j++){
      int i_beg = i - 1;
      int i_end = i + window_size;
      tmp[j] -= geno.at(j, i_beg) * loadings[i_beg] * nSNP/ (window_size * sigma);
      tmp[j] += geno.at(j, i_end) * loadings[i_end] * nSNP/ (window_size * sigma);
      if (direction == 1){
        if ((tmp[j] - scores[j]) > 0){
          stat[i] += (tmp[j] - scores[j]) * (tmp[j] - scores[j]);   
          n_count++;
        }
      } else if (direction == -1){
        if ((tmp[j] - scores[j]) < 0){
          stat[i] += (tmp[j] - scores[j]) * (tmp[j] - scores[j]);   
          n_count++;
        }
      } else if (direction == 0){
        stat[i] += (tmp[j] - scores[j]);   
        n_count++;
      }
    }
    if (direction == 0){
      stat[i] *= stat[i];
    }
    //stat[i] /= sqrt(nIND + 1 - n_count);
  }
  for (int i = (nSNP - window_size); i < nSNP; i++){
    stat[i] = stat[nSNP - window_size - 1];
  }
  return(stat);
}

