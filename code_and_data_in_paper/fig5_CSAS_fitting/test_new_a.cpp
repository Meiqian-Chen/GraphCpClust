// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::List GridSearch1(const NumericVector x, const NumericVector y,
                  const float beta00, const float beta01,
                  const float beta10, const float beta11, const int n) {
  // 变量声明seq1,seq2用于存放网格点坐标
  vec seq1 = linspace(beta00,beta01,n), seq2 = linspace(beta10,beta11,n);
  vec se(n*n);

  mat par(n*n,6), X;
  const colvec cx = as<colvec>(x);
  const colvec cy = as<colvec>(y);
  std::ostream nullstream(0);
  set_cerr_stream(nullstream);
  int T = x.length() - 1;
  mat ct(T-1, n*n);
  int count=0;
  colvec res, coef;
  se.fill(10000);
  
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if(sign(seq1(i))*sign(seq2(j))<0 && seq1(i)<0){
        continue;
      }
      NumericVector fx = pnorm(seq1(i) + seq2(j)*x);
      colvec Fx = as<colvec>(fx);
      if((accu(pow(Fx,2))/T-pow(accu(Fx)/T,2))>0.01){
        X = join_horiz(ones(T-1),Fx.rows(2,T),cy.rows(1,T-1),cy.rows(0,T-2));
        coef = solve(trans(X)*X,trans(X)*cy.rows(2,T));// fit model y ~ X
        res  = cy.rows(2,T) - X*coef;           // residuals
        ct.col(count) = X*coef;
        se(count) = accu(pow(res,2));
        par(count,0) = seq1(i);
        par(count,1) = seq2(j);
        par(count,2) = coef(0);
        par(count,3) = coef(1);
        par(count,4) = coef(2);
        par(count,5) = coef(3);
        count=count+1;
      }
    }
  }
  const uword i = se.index_min();
  return List::create(Named("beta1") = par(i,0),
                      Named("beta2") = par(i,1),
                      Named("a") = par(i,2),
                      Named("b") = par(i,3),
                      Named("alpha1") = par(i,4),
                      Named("alpha2") = par(i,5),
                      Named("residual") = se(i),
                      Named("pred") = ct.col(i));
}
