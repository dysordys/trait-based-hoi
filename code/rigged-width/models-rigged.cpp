// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;


// [[Rcpp::export]]
double errorf(double x) {
  return(2.0*R::pnorm(x*M_SQRT2, 0.0, 1.0, 1, 0) - 1.0);
}


// [[Rcpp::export]]
double sigmoid(double x) {
  return(R::pnorm(x*M_SQRT2, 0.0, 1.0, 1, 0));
}


// [[Rcpp::export]]
double smoothstep(double x) {
  double y;
  if (x > 1.0) y = 1.0; else if (x < 0.0) y = 0.0; else y = x*x*x*(10.0+x*(-15.0+6.0*x));
  return(y);
}


// [[Rcpp::export]]
NumericVector b_rectangular(double theta, NumericVector m) {
  int S = m.size();
  NumericVector b(S);
  for (int i = 0; i < S; i++) b(i) = (m(i) >= -theta) * (m(i) <= theta);
  return(b);
}


// [[Rcpp::export]]
NumericVector b_hier(NumericVector m) {
  int S = m.size();
  NumericVector b(S);
  for (int i = 0; i < S; i++) b(i) = 1.0 - exp(-m(i));
  return(b);
}


// [[Rcpp::export]]
List evo(double time, NumericVector n, List pars) {
  int S = pars["S"], i, j;
  double w = pars["w"], w2 = w*w, alpha_n;
  NumericVector dndt(S), m = pars["m"];
  NumericMatrix alpha(S,S);
  NumericVector b = b_rectangular(pars["theta"], m);
  for (i = 0; i < S; i++) alpha(i,i) = 1.0;
  for (i = 0; i < S - 1; i++) {
    for (j = i + 1; j < S; j++) {
      alpha(i,j) = exp(-(m(i) - m(j)) * (m(i) - m(j)) / w2);
      alpha(j,i) = alpha(i,j);
    }
  }
  for (i = 0; i < S; i++) {
    alpha_n = 0.0;
    for (j = 0; j < S; j++) {
      alpha_n = alpha_n + alpha(i,j)*n(j);
    }
    dndt(i) = n(i) * (b(i) - alpha_n) * smoothstep(n(i)*1.0e6);
  }
  return(List::create(dndt));
}


// [[Rcpp::export]]
List hier(double time, NumericVector n, List pars) {
  int S = pars["S"], i, j;
  double w = pars["w"], alpha_n;
  NumericVector dndt(S), m = pars["m"];
  NumericVector b = b_hier(m);
  for (i = 0; i < S; i++) {
    alpha_n = 0.0;
    for (j = 0; j < S; j++) {
      alpha_n = alpha_n + sigmoid((m(i) - m(j)) / w) * n(j);
    }
    dndt(i) = n(i) * (b(i) - alpha_n) * smoothstep(n(i)*1.0e6);
  }
  return(List::create(dndt));
}


// [[Rcpp::export]]
List evoHOI(double time, NumericVector n, List pars) {
  int S = pars["S"], i, j, k;
  double w = pars["w"], w2 = w*w, kappa = pars["kappa"];
  double dm, dm_ij, dm_jk, dm_ki, num, denom, epsilon_ijk;
  NumericVector dndt(S), m = pars["m"];
  NumericMatrix alpha(S,S);
  NumericVector b = b_rectangular(pars["theta"], m);
  NumericVector alpha_n(S), epsilon_nn(S);
  for (i = 0; i < S; i++) {
    alpha_n(i) = 0.0;
    epsilon_nn(i) = 0.0;
    for (j = 0; j < S; j++) {
      dm = m(i) - m(j);
      alpha(i,j) = exp(-dm*dm/w2);
      alpha_n(i) = alpha_n(i) + alpha(i,j)*n(j);
      for (k = 0; k < S; k++) {
        dm_ij = dm*dm;
        dm_jk = (m(j) - m(k))*(m(j) - m(k));
        dm_ki = (m(k) - m(i))*(m(k) - m(i));
        num = 2.0*w2*(dm_ij + dm_jk + dm_ki);
        denom = 3.0*w2*w2;
        epsilon_ijk = kappa * 2.0*w2*exp(-num/denom) / sqrt(M_SQRT_PI*w*denom);
        epsilon_nn(i) = epsilon_nn(i) + epsilon_ijk*n(j)*n(k);
      }
    }
  }
  for (i = 0; i < S; i++) {
    dndt(i) = n(i) * (b(i) - alpha_n(i) - epsilon_nn(i)) * smoothstep(n(i)*1.0e6);
  }
  return(List::create(dndt));
}


// [[Rcpp::export]]
List evoHOI_rigged(double time, NumericVector n, List pars) {
  int S = pars["S"], i, j, k;
  double w = pars["w"], w2 = w*w, W = pars["W"], W2 = W*W, kappa = pars["kappa"];
  double dm, dm_ij, dm_jk, dm_ki, num, epsilon_ijk;
  NumericVector dndt(S), m = pars["m"];
  NumericMatrix alpha(S,S);
  NumericVector b = b_rectangular(pars["theta"], m);
  NumericVector alpha_n(S), epsilon_nn(S);
  for (i = 0; i < S; i++) {
    alpha_n(i) = 0.0;
    epsilon_nn(i) = 0.0;
    for (j = 0; j < S; j++) {
      dm = m(i) - m(j);
      alpha(i,j) = exp(-dm*dm/w2);
      alpha_n(i) = alpha_n(i) + alpha(i,j)*n(j);
      for (k = 0; k < S; k++) {
        dm_ij = dm*dm;
        dm_jk = (m(j) - m(k))*(m(j) - m(k));
        dm_ki = (m(k) - m(i))*(m(k) - m(i));
        num = 2.0*w2*(dm_ij + dm_jk + dm_ki);
        epsilon_ijk = kappa*exp(-num/W2);
        epsilon_nn(i) = epsilon_nn(i) + epsilon_ijk*n(j)*n(k);
      }
    }
  }
  for (i = 0; i < S; i++) {
    dndt(i) = n(i) * (b(i) - alpha_n(i) - epsilon_nn(i)) * smoothstep(n(i)*1.0e6);
  }
  return(List::create(dndt));
}


// [[Rcpp::export]]
List hierHOI(double time, NumericVector n, List pars) {
  int S = pars["S"], i, j, k;
  double w = pars["w"], W = pars["W"], kappa = pars["kappa"];
  double z0 = pars["z0"], dm, ddm, epsilon_ijk;
  NumericVector dndt(S), m = pars["m"];
  NumericMatrix alpha(S,S);
  NumericVector b = b_hier(m);
  NumericVector alpha_n(S), epsilon_nn(S);
  for (i = 0; i < S; i++) {
    alpha_n(i) = 0.0;
    epsilon_nn(i) = 0.0;
    for (j = 0; j < S; j++) {
      dm = m(i) - m(j);
      alpha(i,j) = sigmoid(dm / w);
      alpha_n(i) = alpha_n(i) + alpha(i,j)*n(j);
      for (k = 0; k < S; k++) {
        ddm = z0 + m(i) - m(j)/2.0 - m(k)/2.0;
        epsilon_ijk = kappa * sigmoid(ddm / W);
        epsilon_nn(i) = epsilon_nn(i) + epsilon_ijk*n(j)*n(k);
      }
    }
  }
  for (i = 0; i < S; i++) {
    dndt(i) = n(i) * (b(i) - alpha_n(i) - epsilon_nn(i)) * smoothstep(n(i)*1.0e6);
  }
  return(List::create(dndt));
}


// [[Rcpp::export]]
List hierHOI_rigged(double time, NumericVector n, List pars) {
  int S = pars["S"], i, j, k;
  double w = pars["w"], W = pars["W"], kappa = pars["kappa"];
  double z0 = pars["z0"], ddm, epsilon_ijk;
  NumericVector dndt(S), m = pars["m"];
  NumericMatrix alpha(S,S);
  NumericVector b = b_hier(m);
  NumericVector alpha_n(S), epsilon_nn(S);
  for (i = 0; i < S; i++) {
    alpha_n(i) = 0.0;
    epsilon_nn(i) = 0.0;
    for (j = 0; j < S; j++) {
      alpha(i,j) = sigmoid((m(i) - m(j)) / w);
      alpha_n(i) = alpha_n(i) + alpha(i,j)*n(j);
      for (k = 0; k < S; k++) {
        ddm = z0 + m(i) - m(j)/2.0 - m(k)/2.0;
        epsilon_ijk = kappa * sigmoid(ddm / W);
        epsilon_nn(i) = epsilon_nn(i) + epsilon_ijk*n(j)*n(k);
      }
    }
  }
  for (i = 0; i < S; i++) {
    dndt(i) = n(i) * (b(i) - alpha_n(i) - epsilon_nn(i)) * smoothstep(n(i)*1.0e6);
  }
  return(List::create(dndt));
}
