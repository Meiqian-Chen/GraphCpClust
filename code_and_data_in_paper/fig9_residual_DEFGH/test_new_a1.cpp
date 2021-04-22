// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::List GridSearch2(const NumericVector x, const arma::colvec& y, 
                  const float beta00, const float beta01,
                  const float beta10, const float beta11, const int n) {
  vec seq1 = linspace(beta00,beta01,n);
  vec seq2 = linspace(beta10,beta11,n);
  vec coef0(n), coef1(n), se(n*n);
  mat par(n*n,4), X;
  std::ostream nullstream(0);
  set_cerr_stream(nullstream);
  int T = x.length();
  int count=0;
  colvec res, coef;
  se.fill(10000);
  
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      NumericVector fx = pnorm(seq1(i) + seq2(j)*x);
      if(sign(seq1(i))*sign(seq2(j))<0 && seq1(i)<0){
        continue;
      }
      colvec Fx = as<colvec>(fx);
      if((accu(pow(Fx,2))/T-pow(accu(Fx)/T,2))>0.01){
        X = join_horiz(ones(size(y)),Fx);
        coef = solve(trans(X)*X,trans(X)*y);// fit model y ~ X
        res  = y - X*coef;           // residuals
        se(count) = accu(pow(res,2));
        par(count,0) = seq1(i);
        par(count,1) = seq2(j);
        par(count,2) = coef(0);
        par(count,3) = coef(1);
        count=count+1;
      }
    }
  }
  const uword i = se.index_min();
  return List::create(Named("beta1") = par(i,0),
                      Named("beta2") = par(i,1),
                      Named("a") = par(i,2),
                      Named("b") = par(i,3),
                      Named("residual") = se(i));
}
