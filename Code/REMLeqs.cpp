#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector REMLeqs(double x, double y, double z, NumericVector R, NumericMatrix Bqq, int M){
  NumericVector out (3);
  for(int i=0; i<M; ++i){
    out[0]+= 2*pow(x,2)*pow(y,2)*(1-pow(z,2)) +z*x*y*(Bqq(i,M+i)+Bqq(M+i,i)+2*R[i]*R[M+i]) - 2*pow(y,2)*(Bqq(i,i)+pow(R[i],2));
    out[1]+= 2*pow(x,2)*pow(y,2)*(1-pow(z,2)) +z*x*y*(Bqq(i,M+i)+Bqq(M+i,i)+2*R[i]*R[M+i]) - 2*pow(x,2)*(Bqq(M+i,M+i)+pow(R[M+i],2));
    out[2]+= -2*z*(1-pow(z,2))*pow(x,2)*pow(y,2) + 2*z*pow(y,2)*(Bqq(i,i)+pow(R[i],2)) + 2*z*pow(x,2)*(Bqq(M+i,M+i)+pow(R[M+i],2)) - 2*x*y*(1+pow(z,2))*(Bqq(i,M+i)+Bqq(M+i,i)+2*R[i]*R[M+i]);
  }
  return out;
}