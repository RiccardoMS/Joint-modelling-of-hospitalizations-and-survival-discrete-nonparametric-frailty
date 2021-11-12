#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector REMLupdate(NumericVector R, NumericMatrix Bqq, int M) {
  NumericVector T(3);
  NumericVector out(3);
  
  for(int i=0;i<M;i++){
    T[0]+=Bqq(i,i)+pow(R[i],2);
    T[2]+=Bqq(M+i,M+i)+pow(R[M+i],2);
    T[1]+=(Bqq(i,M+i)+Bqq(M+i,i)+2*R[i]*R[M+i])/2;
  }
  
  out[1]=T[2]/M;
  out[0]=T[0]/M;
  out[2]=T[1]/sqrt(T[0]*T[2]);
  
  return out;
  
  
}


