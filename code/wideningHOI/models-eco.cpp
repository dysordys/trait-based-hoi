// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;


// [[Rcpp::export]]
double errorf(double x) {
  return 2.0*R::pnorm(x*M_SQRT2, 0.0, 1.0, 1, 0) - 1.0;
}


// [[Rcpp::export]]
double sigmoid(double x) {
  return R::pnorm(x*M_SQRT2, 0.0, 1.0, 1, 0);
}


// [[Rcpp::export]]
double smoothstep(double x) {
  double y;
  if (x > 1.0) y = 1.0; else if (x < 0.0) y = 0.0; else y = x*x*x*(10.0+x*(-15.0+6.0*x));
  return y;
}


// [[Rcpp::export]]
NumericVector b_rectangular(NumericVector s, double theta, NumericVector m) {
  int S = m.size();
  NumericVector b(S);
  for (int i = 0; i < S; i++) {
    b(i) = (errorf((theta+m(i))/(M_SQRT2*s(i)))+errorf((theta-m(i))/(M_SQRT2*s(i))))/2.0;
  }
  return b;
}


// [[Rcpp::export]]
NumericVector b_hier(NumericVector s, double theta, NumericVector m) {
  int S = m.size();
  NumericVector b(S);
  for (int i = 0; i < S; i++) {
    b(i) = 1 - exp(theta + s(i)*s(i)/2 - m(i));
  }
  return b;
}


// [[Rcpp::export]]
List evoHOI(double time, NumericVector state, List pars) {
  int S = pars["S"], i, j, k;
  double w = pars["w"], w2 = w*w, kappa = pars["kappa"], W = pars["W"], W2 = W*W;
  double dm_ij, dm_jk, dm_ki, num, denom, epsilon_ijk;
  NumericVector dndt(S), n = state, m = pars["m"];
  NumericVector sig = pars["sigma"], h2 = pars["h2"], var(S);
  NumericVector b = b_rectangular(sig, pars["theta"], m);
  NumericVector epsilon_nn(S);
  for (i = 0; i < S; i++) var(i) = sig(i)*sig(i);
  for (i = 0; i < S; i++) {
    epsilon_nn(i) = 0.0;
    for (j = 0; j < S; j++) {
      for (k = 0; k < S; k++) {
        dm_ij = (m(i) - m(j))*(m(i) - m(j));
        dm_jk = (m(j) - m(k))*(m(j) - m(k));
        dm_ki = (m(k) - m(i))*(m(k) - m(i));
        num = 8.0*(dm_jk*var(i) + dm_ki*var(j) + dm_ij*var(k)) +
          2.0*w2*(dm_ij + dm_jk + dm_ki);
        denom = W2 + 16.0*(var(i)*var(j) + var(j)*var(k) + var(k)*var(i)) +
          8.0*w2*(var(i) + var(j) + var(k)) + 3.0*w2*w2;
        epsilon_ijk = kappa * 2.0*w2*exp(-num/denom) / sqrt(M_SQRT_PI*w*denom);
        epsilon_nn(i) = epsilon_nn(i) + epsilon_ijk*n(j)*n(k);
      }
    }
  }
  for (i = 0; i < S; i++) {
    dndt(i) = n(i) * (b(i) - epsilon_nn(i)) * smoothstep(n(i)*1.0e6);
  }
  return(List::create(dndt));
}


// [[Rcpp::export]]
List hierHOI(double time, NumericVector state, List pars) {
  int S = pars["S"], i, j, k;
  double W = pars["W"], W2 = W*W, kappa = pars["kappa"];
  double z0 = pars["z0"], svH, ddm, epsilon_ijk;
  NumericVector dndt(S), n = state, m = pars["m"];
  NumericVector sig = pars["sigma"], h2 = pars["h2"], var(S);
  NumericVector b = b_hier(sig, pars["theta"], m);
  NumericVector epsilon_nn(S);
  for (i = 0; i < S; i++) var(i) = sig(i)*sig(i);
  for (i = 0; i < S; i++) {
    epsilon_nn(i) = 0.0;
    for (j = 0; j < S; j++) {
      for (k = 0; k < S; k++) {
        svH = W2 + 2.0*var(i) + var(j)/2.0 + var(k)/2.0;
        ddm = z0 + m(i) - m(j)/2.0 - m(k)/2.0;
        epsilon_ijk = kappa * sigmoid(ddm / sqrt(svH));
        epsilon_nn(i) = epsilon_nn(i) + epsilon_ijk*n(j)*n(k);
      }
    }
  }
  for (i = 0; i < S; i++) {
    dndt(i) = n(i) * (b(i) - epsilon_nn(i)) * smoothstep(n(i)*1.0e6);
  }
  return(List::create(dndt));
}
