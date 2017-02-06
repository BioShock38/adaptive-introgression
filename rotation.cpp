#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

#define NA 9

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat fJ_cpp(int n){
  arma::mat zz(n, n);
  zz.ones();
  zz /= n;
  arma::mat H(n, n);
  H.eye();
  H -= zz;
  return(H);
}

// [[Rcpp::export]]
arma::mat fcnt_cpp(arma::mat &a){
  int nrow = a.n_rows;
  arma::mat aa = fJ_cpp(nrow);
  return(aa * a);
}

// [[Rcpp::export]]
double proc_scale(arma::vec &d, arma::mat &b){
  int K = d.n_elem;
  double scale = 0;
  double ctrace = norm(b, 2);
  for (int k = 0; k < K; k++){
    scale += d[k] / ctrace;
  }
  return(scale);
}

// [[Rcpp::export]]
Rcpp::List pca_rotation(arma::mat &a, arma::mat &b){
  arma::mat fcnt_a = fcnt_cpp(a);
  arma::mat fcnt_b = fcnt_cpp(b);
  arma::mat x = fcnt_a.t() * fcnt_b;
  arma::mat u;
  arma::vec s;
  arma::mat v;
  svd(u, s, v, x);
  double scale = proc_scale(s, b);
  arma::mat R;
  R = v * u.t();
  arma::cx_vec eigval_v;
  arma::cx_mat eigvec_v;
  arma::cx_vec eigval_u;
  arma::cx_mat eigvec_u;
  eig_gen(eigval_v, eigvec_v, v);
  eig_gen(eigval_u, eigvec_u, u);
  arma::cx_double tmp_v = prod(eigval_v);
  arma::cx_double tmp_u = prod(eigval_u);
  double chk1 = real(tmp_v);
  double chk2 = real(tmp_u);
  if ((chk1 < 0) && (chk2 > 0)) {
    for (int i = 0; i < v.n_rows; i++){
      v(i, v.n_cols - 1) *= (-1);
    }
    R = v * u.t();
  }
  if ((chk2 < 0) && (chk1 > 0)) {
    for (int i = 0; i < u.n_rows; i++){
      u(i, u.n_cols - 1) *= (-1);
    }
    R = v * u.t();
  }
  return Rcpp::List::create(Rcpp::Named("R") = R,
                            Rcpp::Named("s") = scale);
}

// [[Rcpp::export]]
arma::vec cmpt_local_scores_0(arma::mat &x, arma::mat &u, arma::mat &v, arma::vec &s, int window_size, int direction){
  int nIND = x.n_rows;  
  int K = u.n_cols;
  arma::vec sum_loadings_sq(K);
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      for (int p = 0; p < window_size; p++){
        u(j, k) += x.at(j, p) * v(p, k) / s[k];
        sum_loadings_sq[k] += v(p, k) * v(p, k);
      }      
    }
  }
  return(sum_loadings_sq);
}

// [[Rcpp::export]]
arma::mat rescale_local_scores(arma::mat &uloc, arma::vec &slsq){
  arma::mat uscaled(uloc.n_rows, uloc.n_cols);
  for (int j = 0; j < uloc.n_rows; j++){
    for (int k = 0; k < uloc.n_cols; k++){
      uscaled(j, k) = uloc(j, k) / slsq[k];
    }
  }
  return(uscaled);
}

// [[Rcpp::export]]
double cmpt_proc_dist(arma::mat &xref, arma::mat &xloc){
  arma::mat resid = xref - xloc;
  arma::mat txx = resid.t() * resid;
  return(sum(txx.diag()));
}

// [[Rcpp::export]]
double cmpt_dist_i(arma::mat &uref, arma::vec &uloc, int K, int direction){
  double res = 0;
  for (int j = 0; j < uref.n_rows; j++){
    //double tmp = uloc(j, K) - uref(j, K); 
    double tmp = uloc[j] - uref(j, K); 
    if (direction == 1){
      if (tmp > 0){
        res += tmp * tmp;
      }
    } else if (direction == (-1)){
      if (tmp < 0){
        res += tmp * tmp;
      }  
    } else if (direction == 0){
      res += tmp; 
    }
  }
  return(res);
}

// [[Rcpp::export]]
void cmpt_local_scores_i(arma::mat &x, arma::mat &u, arma::mat &v, arma::vec &s, int window_size, int direction, int i, arma::vec &slsq){
  int nIND = x.n_rows;  
  int K = u.n_cols;
  int i_beg = i - 1;
  int i_end = i + window_size; 
  for (int j = 0; j < nIND; j++){
    for (int k = 0; k < K; k++){
      u(j, k) -= (x.at(j, i_beg) * v(i_beg, k)) / s[k];
      u(j, k) += (x.at(j, i_end) * v(i_end, k)) / s[k];
    }
  }
  for (int k = 0; k < K; k++){
    slsq[k] -= (v(i_beg, k) * v(i_beg, k));
    slsq[k] += (v(i_end, k) * v(i_end, k));
  }
}

// [[Rcpp::export]]
arma::vec cmpt_ms(arma::mat &u_anc_global, arma::mat &u_anc_local, int axis){
  double m_global = mean(u_anc_global.col(axis));
  double m_local = mean(u_anc_local.col(axis));
  double s_global = stddev(u_anc_global.col(axis));
  double s_local = stddev(u_anc_local.col(axis));
  arma::vec res(2);
  res[0] = m_global - m_local;
  res[1] = s_global / s_local;
  return(res);
}

// [[Rcpp::export]]
arma::vec align_scores(arma::mat &u_anc_global, arma::mat &u_anc_local, arma::mat &u_adm_local, int axis){
  int N = u_adm_local.n_rows;
  arma::vec res = cmpt_ms(u_anc_global, u_anc_local, axis);
  arma::vec us(N);
  for (int n = 0; n < N; n++){
    us[n] = (u_adm_local(n, axis) + res[0]) * res[1];
    //us[n] = (u_adm_local(n, axis));
  }
  return(us);
}


// [[Rcpp::export]]
arma::vec compute_stat_5(arma::mat &geno, 
                          arma::mat &scores, 
                          arma::mat &anc_geno, 
                          arma::mat &anc_scores, 
                          int PC, 
                          arma::mat &loadings, 
                          arma::vec sigma, 
                          int window_size, 
                          int direction){
  int nSNP = geno.n_cols;
  int nIND = geno.n_rows;
  int nANC = anc_geno.n_rows;
  int K = scores.n_cols;
  int axis = PC - 1;
  arma::vec stat(nSNP);
  stat.zeros();
  arma::mat uloc(nIND, K);
  uloc.zeros();
  arma::mat ancloc(nANC, K);
  ancloc.zeros();
  arma::vec slsq(K);
  arma::vec slsq2(K);
  arma::vec uscaled(nIND);
  
  slsq = cmpt_local_scores_0(geno, uloc, loadings, sigma, window_size, direction);
  slsq2 = cmpt_local_scores_0(anc_geno, ancloc, loadings, sigma, window_size, direction);
  uscaled = align_scores(anc_geno, ancloc, uloc, axis);
  //uscaled *= nSNP / window_size;
  //uscaled /= slsq[0];
  stat[0] = cmpt_dist_i(scores, uscaled, axis, direction);
  
  for (int i = 1; i < (nSNP - window_size); i++){
    cmpt_local_scores_i(geno, uloc, loadings, sigma, window_size, direction, i, slsq);
    cmpt_local_scores_i(anc_geno, ancloc, loadings, sigma, window_size, direction, i, slsq);
    uscaled = align_scores(anc_geno, ancloc, uloc, axis);
    //uscaled *= nSNP / window_size;
    //uscaled /= slsq[0];
    stat[i] = cmpt_dist_i(scores, uscaled, axis, direction);
  }
  for (int i = (nSNP - window_size); i < nSNP; i++){
    stat[i] = stat[nSNP - window_size - 1];
  }
  return(stat);
}

 
// // [[Rcpp::export]]
// Rcpp::List compute_stat_3(arma::mat &geno, 
//                           arma::mat &scores, 
//                           arma::mat &anc_geno, 
//                           arma::mat &anc_scores, 
//                           int PC, 
//                           arma::mat &loadings, 
//                           arma::vec sigma, 
//                           int window_size, 
//                           int direction){
//   int nSNP = geno.n_cols;
//   int nIND = geno.n_rows;
//   int nANC = anc_geno.n_rows;
//   int K = scores.n_cols;
//   int axis = PC - 1;
//   arma::vec stat(nSNP);
//   stat.zeros();
//   arma::vec stat_proc(nSNP);
//   stat_proc.zeros();
//   arma::vec tmp(nIND);
//   tmp.zeros();
//   arma::mat uloc(nIND, K);
//   uloc.zeros();
//   arma::mat ancloc(nANC, K);
//   ancloc.zeros();
//   arma::vec slsq(K);
//   arma::mat uscaled(nIND, K);
//   arma::mat ancscaled(nIND, K);
//   
//   
//   slsq = cmpt_local_scores_0(geno, uloc, loadings, sigma, window_size, direction);
//   slsq = cmpt_local_scores_0(anc_geno, ancloc, loadings, sigma, window_size, direction);
//   // for (int k = 0; k < K; k++){
//   //   slsq[k] = double(window_size) / double(nSNP);
//   // }
//   uscaled = rescale_local_scores(uloc, slsq);
//   ancscaled = rescale_local_scores(ancloc, slsq);
//   stat[0] = cmpt_dist_i(scores, uscaled, axis, direction);
// 
//   for (int i = 1; i < (nSNP - window_size); i++){
//     cmpt_local_scores_i(geno, uloc, loadings, sigma, window_size, direction, i, slsq);
//     cmpt_local_scores_i(anc_geno, ancloc, loadings, sigma, window_size, direction, i, slsq);
//     // for (int k = 0; k < K; k++){
//     //   slsq[k] = double(window_size) / double(nSNP);
//     // }
//     uscaled = rescale_local_scores(uloc, slsq);
//     ancscaled = rescale_local_scores(ancloc, slsq);
//     arma::mat tt = pca_rotation(anc_scores, ancscaled);
//     arma::mat Bhat = uscaled * tt;
//     stat[i] = cmpt_dist_i(scores, uscaled, axis, direction);
//     stat_proc[i] = cmpt_dist_i(scores, Bhat, axis, direction);
//   }
//   for (int i = (nSNP - window_size); i < nSNP; i++){
//     stat[i] = stat[nSNP - window_size - 1];
//     stat_proc[i] = stat_proc[nSNP - window_size - 1];
//   }
//   return Rcpp::List::create(Rcpp::Named("stat") = stat,
//                             Rcpp::Named("stat_proc") = stat_proc);
// }

// [[Rcpp::export]]
arma::mat match_angle(arma::vec A, arma::vec B, arma::vec C, arma::vec D){
  arma::vec AB(2);
  arma::vec CD(2);
  AB[0] = B[0] - A[0];
  AB[1] = B[1] - A[1];
  CD[0] = D[0] - C[0];
  CD[1] = D[1] - C[1];
  double ABdotCD = arma::dot(AB, CD);
  double normAB = arma::norm(AB, 2);
  double normCD = arma::norm(CD, 2);
  double theta = acos(ABdotCD / (normAB * normCD));
  arma::mat R(2, 2);
  R(0, 0) = cos(theta);
  R(0, 1) = (-1) * sin(theta);
  R(1, 0) = sin(theta);
  R(1, 1) = cos(theta);
  Rprintf("%f\n", theta);
  return(R);
}

// [[Rcpp::export]]
double match_length(arma::vec A, arma::vec B, arma::vec C, arma::vec D){
  arma::vec AB(2);
  arma::vec CD(2);
  AB[0] = B[0] - A[0];
  AB[1] = B[1] - A[1];
  CD[0] = D[0] - C[0];
  CD[1] = D[1] - C[1];
  double normAB = arma::norm(AB, 2);
  double normCD = arma::norm(CD, 2);
  double s = normCD / normAB;
  return(s);
}

// [[Rcpp::export]]
arma::vec cmpt_centroid(arma::mat &x){
  int dim = x.n_cols;
  int n_obs = x.n_rows;
  arma::vec a(dim);
  a.zeros();
  for (int d = 0; d < dim; d++){
    for (int i = 0; i < n_obs; i++){
      a[d] += double (x(i, d) / n_obs);  
    }
  }
  return(a);
}
