#include <RcppArmadillo.h>
#include <Rcpp.h>
//#include "/home/shushuz/R_libs/RcppArmadillo/include/RcppArmadillo.h"
//#include "/home/shushuz/R_libs/Rcpp/include/Rcpp.h"
#include <math.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec rho_rq_function(vec u, double theta) {
    return theta * arma::max(u, vec(u.n_elem,fill::zeros)) + (theta - 1) * arma::min(u, vec(u.n_elem,fill::zeros));
}
// [[Rcpp::export]]
vec quickrq(const vec& y, const mat& x, const double& tau=0.5, const int& kmax=20, const double& big=1e+20, const double& eps=0.000001, const double& beta=0.97) {
    int n = y.n_elem;
    int p = x.n_cols;
    //if(int_ == true) {
    //    mat ones(n, 1, fill::ones);
    //    x.insert_cols(0, ones);
    //}
    //if(w.n_elem == 0) {
    vec w = vec(n, fill::zeros);
    //}

    double yw = big;
    int k = 1;
    vec coef(p);
    while(k <= kmax && std::abs(dot(y, w)/yw - 1) > eps) {
        vec d = arma::min(tau - w, 1 - tau + w);
        //arma::mat W = arma::diagmat(d);
	//arma::mat X_d = x.each_col() % d;
        //vec y_d = y % d;
  // Solve the weighted least squares problem
	coef = arma::solve(x.each_col() % d, y % d);
//Rcout<<"k:"<<k<<endl;
//	Rcout<<coef<<endl;
  	vec resid = y - x * coef;
        yw = sum(vec(rho_rq_function(resid, tau)));
        k++;
        vec s =  resid % d % d;
	double alpha = std::max(eps, std::max(static_cast<double>(arma::max(s/(tau - w))), static_cast<double>(arma::max(-s/(1 - tau + w)))));
	//double alpha = std::max(arma::datum::eps, arma::max(arma::max(s/(tau - w)),arma::max( -s/(1 - tau + w))))[0];
        w = w + (beta/alpha) * s;
    }
 //   Rcout<<"kmax:"<<k<<endl;
    //vec res = y - x * beta;
    //double zero_res = arma::sort(arma::abs(res))(X.n_rows - 1); 
    //uvec zero_id = arma::find(abs(res) <= zero_res);
    //for(arma::uword j=0; j<zero_id.n_elem; j++) {
    //    res(zero_id(j)) = 0;
    //}
    return coef;
    //return List::create(Named("coef")=coef, Named("res")=res);
}


// [[Rcpp::export]]
int Match_Strings(StringVector S, String s) {
  for(int i = 0; i< S.length();i++){
    if(S(i) == s) return(i);
  }
  return -1;
}
// [[Rcpp::export]]
uvec my_interpolate(vec tau_all, vec tau_use){
  // From a really find grid tau_all, find the index of tau_use
  int T = tau_use.n_elem;
  uvec idx_target(T);
    int idx_grid = 0;
    double eps = 1;
    for(int t = 0; t < T; t++){
      while((idx_grid < tau_all.n_elem) & (tau_all(idx_grid) < tau_use(t)) ){
        eps = tau_use(t) - tau_all(idx_grid);
        idx_grid ++;
      }
      if (idx_grid < tau_all.n_elem){
        if(eps > tau_all(idx_grid)- tau_use(t)){
          idx_target(t) = idx_grid;
        } else {
          idx_target(t) = idx_grid - 1;
        }
        
      } else {
        // Reached the end...
        idx_target(t,T-1) = tau_all.n_elem-1;
        break;
      }
 } 
 return idx_target;
}



// [[Rcpp::export]]
rowvec my_lm(mat X, mat y, rowvec pred_x) {
  /* Function to compute linear regression
   * Input X: n by (p+1), includes an intercept term
   * Input Y: n by q, can be multi-dimensional
   * Input pred_x: 1 by (p+1), includes an intercept term
   * Output: a row vector 1 by q; the predicted values at pred_x
  */
  rowvec a;
  mat SX = X.t()* X;
  if(rcond(SX) < 1e-10){ // if design is almost singular
    a = mean(y,0);
  } else {
    mat beta = solve(X.t()* X, X.t()*y);
    a = pred_x*beta;
  }
  return a;
}

// Gaussian Kernel function
double gaussian_kernel(double x, double mu, double bandwidth) {
  return exp(-0.5 * pow((x - mu) / bandwidth, 2)) / (bandwidth * sqrt(2 * M_PI));
}

vec kde(vec data, vec evaluation_points, double bandwidth) {
 vec result(evaluation_points.n_elem, fill::zeros);
  
  for (size_t i = 0; i < evaluation_points.n_elem; ++i) {
    for (size_t j = 0; j < data.n_elem; ++j) {
      result(i) += gaussian_kernel(evaluation_points(i), data(j), bandwidth);
    }
  }
  
  result /= (data.n_elem);
  
  return result;
}


// [[Rcpp::export]]
mat Local_Qt_KNN(const mat& X, const vec& Y,
                 const vec tau_list){
  /* Function to compute fitted value of (local-) quantile regression at potentially a fine grid of quantile levels
   * Specifically tuned to speed up computation, in comparison to the quantreg package
   * Input X: INCLUDES intercept
   * Input Y: A one-dimensional vector
   * Input tau_list: a vector of targeted quantile levels
   * Output: a matrix of dimension NROW(X) by LENGTH(tau_list) for the fitted values
   */
  int n = X.n_rows;
  int p = X.n_cols;
  int T = tau_list.n_elem;
  mat qt_coef(p,T);
  mat Res_val(n,T);
  vec u(Y.n_elem);
  vec f(3);
  double f0_prime;
  mat Q(p,p);
  vec XXQX(p,fill::zeros);
  mat XXX = zeros(p,p*p);
  vec XX(p*p,fill::zeros);
  vec pt(3,fill::zeros);
  pt[0] = -0.1;
  pt[2] = 0.1;
for(int i=0; i < T; i++){
	//qt_coef.col(i) = fitQuantRegLP(Y,X,tau_list[i]);
	//qt_coef.col(i) = quantreg(Y,X,tau_list[i]);
	qt_coef.col(i) = quickrq(Y,X,tau_list[i]);
	//u = Y - X*qt_coef.col(i);
	//f = kde(u,pt,0.1);
	//f0_prime = (f[2] - f[0])/(pt[2]-pt[0]);
	//Q = arma::inv(X.t()*X/n)/f[1];
	//for (int j = 0; j < n; ++j) {
	//	rowvec X_j = X.row(j);
	//	mat outer_Xj = X_j.t() * X_j;
	//	XXQX +=  outer_Xj * Q * X_j.t();
	//	XXX += arma::kron(outer_Xj, X_j);
	//	XX += arma::kron(X_j.t(), X_j.t());
	//}
	//XXQX /= n;
  	//XXX /= n;
  	//XX /= n;
	//qt_coef.col(i) -= 1/n*Q*((1/2-tau_list[i])*f[1]*XXQX - tau_list[i]*(1-tau_list[i])/2*f0_prime*XXX*arma::kron(Q,Q)*XX);
	//Rcout<<"bias:"<<1/n*Q*((1/2-tau_list[i])*f[1]*XXQX - tau_list[i]*(1-tau_list[i])/2*f0_prime*XXX*arma::kron(Q,Q)*XX)<<endl;
}
Res_val = X*qt_coef;
return Res_val;
}

mat Local_Qt_KNN_coef(const mat& X, const vec& Y,
                 const vec tau_list){
  /* Function to compute fitted value of (local-) quantile regression at potentially a fine grid of quantile levels
   * Specifically tuned to speed up computation, in comparison to the quantreg package
   * Input X: INCLUDES intercept
   * Input Y: A one-dimensional vector
   * Input tau_list: a vector of targeted quantile levels
   * Output: a matrix of dimension NROW(X) by LENGTH(tau_list) for the fitted values
   */
  int n = X.n_rows;
  int p = X.n_cols;
  int T = tau_list.n_elem;
  mat qt_coef(p,T);

for(int i=0; i < T; i++){
        //qt_coef.col(i) = fitQuantRegLP(Y,X,tau_list[i]);
        //qt_coef.col(i) = quantreg(Y,X,tau_list[i]);
        qt_coef.col(i) = quickrq(Y,X,tau_list[i]);
}
return qt_coef;
}

// [[Rcpp::export]]
rowvec Local_SQ_Neyman(const mat& X, const vec& Y, const mat& Qt, 
                     const vec& tau_grid, const rowvec& pred_x){
  /* Function to evaluate the predicted value at one point for the 
   * (local-)linear SQ estimation using the Neyman-Orthogonalized score function
   * Input X: a matrix that INCLUDES intercept
   * Input Y: a one-dimensional vector
   * Input Qt: pre-estimated conditional quantile matrix of dimension NROW(X) by LENGTH(tau_grid)
   * Input tau_grid: a vector of targeting quantile levels
   * Input pred_x: INCLUDES intercept
   * Output: a vector of length 1 by LENGTH(tau_grid); the predicted conditional SQ process at pred_x
   */
  int K = Y.n_elem;
  int T = tau_grid.n_elem;
  mat Zlocal(K,T);
  for(int k = 0; k < K; k++){
    for(int t=0;t<T;t++){
      Zlocal(k,t) = (Y(k) - Qt(k,t))*(Y(k) >= Qt(k,t))/(1-tau_grid(t)) + Qt(k,t);
    }
  }
    return(my_lm(X,Zlocal,pred_x));
}

// [[Rcpp::export]]
rowvec Local_SQ_LS(const mat& X, const vec& Y, const mat& Qt, 
                       const vec& tau_grid, const rowvec& pred_x){
  /* Function to evaluate the predicted value at one point for the 
   * (local-)linear SQ estimation using the LS score function
   * Input X: a matrix that INCLUDES intercept
   * Input Y: a one-dimensional vector
   * Input Qt: pre-estimated conditional quantile matrix of dimension NROW(X) by LENGTH(tau_grid)
   * Input tau_grid: a vector of targeting quantile levels
   * Input pred_x: INCLUDES intercept
   * Output: a vector of length 1 by LENGTH(tau_grid); the predicted conditional SQ process at pred_x
   */
  int T = tau_grid.n_elem;
  uvec my_idx;
  rowvec res(T);
  int p = X.n_cols;
  for(int t=0;t<T;t++){
    my_idx = find(Y - Qt.col(t) >= 0); // above tau-th quantile
    if(my_idx.n_elem == 0){ // not enough y above the quantile 
      res(t) = Qt(0,t);
    } else if(my_idx.n_elem <= p){ // not enough y above the quantile 
      res(t) = mean(Y.elem(my_idx));
    } else{
      res(t) = as_scalar(my_lm(X.rows(my_idx), Y.elem(my_idx),pred_x));
    }
  }
  return(res);
}




// [[Rcpp::export]]
mat Rock_KNN0_subsample(const mat& xdata, const vec& ydata, 
              const umat& nbr_idx, const uvec use_idx, vec tau_grid, String sq_score,
              const mat& qt_estimates) {
  /* Main function, computes all required initial SQ estimates based on kNN bins
   * Relies on given quantile regression estimators
   * Input xdata: full data covariate, does NOT contain intercept
   * Input ydata: full data response vector
   * Input use_idx: index vector of subsampled data; length = m
   * Input nbr_idx: m by 'K' matrix, the indices for the K-neighbours of the subsampled data; 
   *                K is number of neighbours
   * Input tau_grid: the grid of quantile levels
   * Input sq_score: one of {'Neyman','LS', 'Both'}
   * Input qt_estimates: the conditional quantile estimates for full data and all quantile levels;
   *                     of dimension NROW(xdata) by LENGTH(tau_grid)
   * Output: a matrix of initial SQ estimates of dimension m by LENGTH(tau_grid);
   *         dimension at 2m by LENGTH(tau_grid) if sq_score = 'Both'
   */
  int n = xdata.n_rows;
  int m  = use_idx.n_elem;
  int p = xdata.n_cols;
  int K = nbr_idx.n_cols;
  int T = tau_grid.n_elem;
  mat X_reg(n, p+1, fill::ones);  // add the intercept
  X_reg.cols(1,p) = xdata;
  uvec v;
  urowvec u(K);
  mat Res(m,T);
  StringVector sq_options = {"LS", "Neyman","Both"};
  int idx_sq = Match_Strings(sq_options,sq_score) + 1;
  
  switch (idx_sq){
  case 2:{ // Neyman SQ
    for(int i = 0; i< m; i++){
    u = nbr_idx.row(i);
    v = u.elem(find(u > 0)) - 1;
    Res.row(i) = Local_SQ_Neyman(X_reg.rows(v),ydata.elem(v),
            qt_estimates.rows(v),tau_grid,X_reg.row(use_idx(i) - 1));
  }
  return(Res);
    break;
  }
  case 1:{ // LS
    for(int i = 0; i< m; i++){
    u = nbr_idx.row(i);
    v = u.elem(find(u > 0)) - 1;
    Res.row(i) = Local_SQ_LS(X_reg.rows(v),ydata.elem(v),
            qt_estimates.rows(v),tau_grid,X_reg.row(use_idx(i) - 1));
  }
  return(Res);
    break;
  }
  case 3:{ // Both
    mat Res2(m,T);
    for(int i = 0; i< m; i++){
    u = nbr_idx.row(i);
    v = u.elem(find(u > 0)) - 1;
    Res.row(i) = Local_SQ_Neyman(X_reg.rows(v),ydata.elem(v),
            qt_estimates.rows(v),tau_grid,X_reg.row(use_idx(i) - 1));
    Res2.row(i) = Local_SQ_LS(X_reg.rows(v),ydata.elem(v),
            qt_estimates.rows(v),tau_grid,X_reg.row(use_idx(i) - 1));
  }
    return(join_cols(Res,Res2));
    break;
  }
  default:{
    stop("Error, incorrect SQ estimation method! \n");
  }
  }
  return Res;
}

// [[Rcpp::export]]
mat quant_res(const mat& xdata, const vec& ydata,
              const umat& nbr_idx,const uvec use_idx, vec tau_grid) {
int n = xdata.n_rows;
  int m = use_idx.n_elem;
  int p = xdata.n_cols;
  int K = nbr_idx.n_cols;
  int T = tau_grid.n_elem;
  mat X_reg(n, p+1, fill::ones);  // include the intercept
  X_reg.cols(1,p) = xdata;
  uvec v;
  urowvec u(K);
  mat Res(m,T), qt_coef(p+1,T);
 // Neyman SQ
    for(int i = 0; i< m; i++){
    u = nbr_idx.row(i);
    v = u.elem(find(u > 0)) - 1;
    qt_coef = Local_Qt_KNN_coef(X_reg.rows(v),ydata.elem(v), tau_grid);
    //Rcout<<"use_idx:"<<use_idx(i)-1<<endl;
    //mat fitted = xdata.rows(v)*qt_coef.rows(1,p);
    Rcout<<"ydata:"<<ydata[use_idx(i)-1]<<endl;
    Rcout<<"fitted:"<<xdata.row(use_idx(i)-1)*qt_coef.rows(1,p)<<endl;
    Res.row(i) = ydata[use_idx(i)-1] - xdata.row(use_idx(i)-1)*qt_coef.rows(1,p);
Rcout<<"Res:"<<Res.row(i)<<endl;
    //Res.row(i) = ydata.elem(v).t() - fitted;
  }
    //Rcout<<"Res:"<<Res<<endl;
    return(Res);
}

// [[Rcpp::export]]
mat Rock_KNN2_subsample(const mat& xdata, const vec& ydata, 
              const umat& nbr_idx,const uvec use_idx, vec tau_grid, String sq_score) {
  /* Main function 2, Same as 'Rock_KNN0_subsample'; but relies on local kNN quantile regression estimators
   * Input xdata: full data covariate, does NOT contain intercept
   * Input ydata: full data response vector
   * Input use_idx: index vector of subsampled data; length = m
   * Input nbr_idx: m by 'K' matrix, the indices for the K-neighbours of the subsampled data; 
   *                K is number of neighbours
   * Input tau_grid: the grid of quantile levels
   * Input sq_score: one of {'Neyman','LS', 'Both'}
   * Output: a matrix of initial SQ estimates of dimension m by LENGTH(tau_grid);
   *         dimension at 2m by LENGTH(tau_grid) if sq_score = 'Both'
   */
  int n = xdata.n_rows;
  int m = use_idx.n_elem;
  int p = xdata.n_cols;
  int K = nbr_idx.n_cols;
  int T = tau_grid.n_elem;
  mat X_reg(n, p+1, fill::ones);  // include the intercept
  X_reg.cols(1,p) = xdata;
  uvec v;
  urowvec u(K);
  mat Res(m,T), qt_bin(K,T);
  StringVector sq_options = {"LS", "Neyman","Both"};
  int idx_sq = Match_Strings(sq_options,sq_score) + 1;
  
  switch (idx_sq){
  case 2:{ // Neyman SQ
    for(int i = 0; i< m; i++){
    u = nbr_idx.row(i);
    v = u.elem(find(u > 0)) - 1;
    qt_bin = Local_Qt_KNN(X_reg.rows(v),ydata.elem(v), tau_grid);
    
    Res.row(i) = Local_SQ_Neyman(X_reg.rows(v),ydata.elem(v),
            qt_bin,tau_grid,X_reg.row(use_idx(i) - 1));
  }
    return(Res);
    break;
  }
  case 1:{ // LS
    for(int i = 0; i< m; i++){
    u = nbr_idx.row(i);
    v = u.elem(find(u > 0)) - 1;
    qt_bin = Local_Qt_KNN(X_reg.rows(v),ydata.elem(v), tau_grid);
    
    Res.row(i) = Local_SQ_LS(X_reg.rows(v),ydata.elem(v),
            qt_bin,tau_grid,X_reg.row(use_idx(i) - 1));
  }
    return(Res);
    break;
  }
  case 3:{ // Return both
    mat Res2(m,T);
    for(int i = 0; i< m; i++){
    u = nbr_idx.row(i);
    v = u.elem(find(u > 0)) - 1;
    qt_bin = Local_Qt_KNN(X_reg.rows(v),ydata.elem(v), tau_grid);
    
    Res.row(i) = Local_SQ_Neyman(X_reg.rows(v),ydata.elem(v),
            qt_bin,tau_grid,X_reg.row(use_idx(i) - 1));
    Res2.row(i) = Local_SQ_LS(X_reg.rows(v),ydata.elem(v),
            qt_bin,tau_grid,X_reg.row(use_idx(i) - 1));
  }
    return(join_cols(Res,Res2));
    break;
  }
  default:{
    stop("Error, incorrect SQ estimation method! \n");
  }
  }
  
  return(Res);
}

