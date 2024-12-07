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
List b_rectangular(NumericVector s, double theta, NumericVector m) {
  int S = m.size();
  double v;
  NumericVector b(S), g(S);
  for (int i = 0; i < S; i++) {
    v = s(i)*s(i);
    b(i) = (errorf((theta+m(i))/(M_SQRT2*s(i))) + errorf((theta-m(i)) /
      (M_SQRT2*s(i))))/2.0;
    g(i) = (exp(-(m(i)+theta)*(m(i)+theta)/(2.0*v)) -
      exp(-(m(i)-theta)*(m(i)-theta)/(2.0*v))) * M_1_SQRT_2PI * s(i);
  }
  return(List::create(b, g));
}


// [[Rcpp::export]]
List b_hier(NumericVector s, double theta, NumericVector m) {
  int S = m.size();
  double v;
  NumericVector b(S), g(S);
  for (int i = 0; i < S; i++) {
    v = s(i)*s(i);
    b(i) = 1 - exp(theta + v/2 - m(i));
    g(i) = v * exp(theta + v/2 - m(i));
  }
  return(List::create(b, g));
}


// [[Rcpp::export]]
List evo(double time, NumericVector state, List pars) {
  int S = pars["S"], i, j;
  double w = pars["w"], w2 = w*w, dm, sv, alpha_n, beta_n, betafact;
  NumericVector dvdt(2*S), n = state[Range(0, S - 1)], m = state[Range(S, 2*S - 1)];
  NumericVector sig = pars["sigma"], h2 = pars["h2"], var(S);
  NumericMatrix alpha(S,S), beta(S,S);
  List bg = b_rectangular(sig, pars["theta"], m);
  NumericVector b = bg[0], g = bg[1];
  for (i = 0; i < S; i++) {
    var(i) = sig(i)*sig(i);
    alpha(i,i) = w/sqrt(4.0*var(i) + w2);
    beta(i,i) = 0.0;
  }
  for (i = 0; i < S - 1; i++) {
    for (j = i + 1; j < S; j++) {
      dm = m(i) - m(j);
      sv = 2.0*(var(i) + var(j));
      alpha(i,j) = exp(-dm*dm/(sv + w2))*w / sqrt(sv + w2);
      alpha(j,i) = alpha(i,j);
      betafact = -alpha(i,j)*2.0*dm / (sv + w2);
      beta(i,j) = betafact*var(i);
      beta(j,i) = -betafact*var(j);
    }
  }
  for (i = 0; i < S; i++) {
    alpha_n = 0.0;
    beta_n = 0.0;
    for (j = 0; j < S; j++) {
      alpha_n = alpha_n + alpha(i,j)*n(j);
      beta_n = beta_n + beta(i,j)*n(j);
    }
    dvdt(i)   = n(i)  * (b(i) - alpha_n) * smoothstep(n(i)*1.0e6);
    dvdt(i+S) = h2(i) * (g(i) -  beta_n) * smoothstep(n(i)*1.0e5);
  }
  return(List::create(dvdt));
}


// [[Rcpp::export]]
List hier(double time, NumericVector state, List pars) {
  int S = pars["S"], i, j;
  double w = pars["w"], w2 = w*w, dm, sv, alpha_n, beta_n;
  NumericVector dvdt(2*S), n = state[Range(0, S - 1)], m = state[Range(S, 2*S - 1)];
  NumericVector sig = pars["sigma"], h2 = pars["h2"];
  List bg = b_hier(sig, pars["theta"], m);
  NumericVector b = bg[0], g = bg[1];
  for (i = 0; i < S; i++) {
    alpha_n = 0.0;
    beta_n = 0.0;
    for (j = 0; j < S; j++) {
      dm = m(i) - m(j);
      sv = 2.0*(sig(i)*sig(i) + sig(j)*sig(j)) + w2;
      alpha_n = alpha_n + sigmoid(dm / sqrt(sv)) * n(j);
      beta_n = beta_n + (sig(i)*sig(i)*exp(-dm*dm / sv) / sqrt(sv*M_PI)) * n(j);
    }
    dvdt(i)   = n(i)  * (b(i) - alpha_n) * smoothstep(n(i)*1.0e6);
    dvdt(i+S) = h2(i) * (g(i) -  beta_n) * smoothstep(n(i)*1.0e5);
  }
  return(List::create(dvdt));
}


// [[Rcpp::export]]
List evoHOI(double time, NumericVector state, List pars) {
  int S = pars["S"], i, j, k;
  double w = pars["w"], w2 = w*w, kappa = pars["kappa"];
  double sv, dm, dm_ij, dm_jk, dm_ki, num, denom, epsilon_ijk, gamma_ijk, numgamma;
  NumericVector dvdt(2*S), n = state[Range(0, S - 1)], m = state[Range(S, 2*S - 1)];
  NumericVector sig = pars["sigma"], h2 = pars["h2"], var(S);
  NumericMatrix alpha(S,S), beta(S,S);
  List bg = b_rectangular(sig, pars["theta"], m);
  NumericVector b = bg[0], g = bg[1];
  NumericVector alpha_n(S), beta_n(S), epsilon_nn(S), gamma_nn(S);
  for (i = 0; i < S; i++) var(i) = sig(i)*sig(i);
  for (i = 0; i < S; i++) {
    alpha_n(i) = 0.0;
    beta_n(i) = 0.0;
    epsilon_nn(i) = 0.0;
    gamma_nn(i) = 0.0;
    for (j = 0; j < S; j++) {
      dm = m(i) - m(j);
      sv = 2.0*(var(i) + var(j));
      alpha(i,j) = exp(-dm*dm/(sv + w2))*w / sqrt(sv + w2);
      beta(i,j) = -alpha(i,j)*2.0*var(i)*dm / (sv + w2);
      alpha_n(i) = alpha_n(i) + alpha(i,j)*n(j);
      beta_n(i) = beta_n(i) + beta(i,j)*n(j);
      for (k = 0; k < S; k++) {
        dm_ij = dm*dm;
        dm_jk = (m(j) - m(k))*(m(j) - m(k));
        dm_ki = (m(k) - m(i))*(m(k) - m(i));
        num = 8.0*(dm_jk*var(i) + dm_ki*var(j) + dm_ij*var(k)) +
          2.0*w2*(dm_ij + dm_jk + dm_ki);
        denom = 16.0*(var(i)*var(j) + var(j)*var(k) + var(k)*var(i)) +
          8.0*w2*(var(i) + var(j) + var(k)) + 3.0*w2*w2;
        epsilon_ijk = kappa * 2.0*w2*exp(-num/denom) / sqrt(M_SQRT_PI*w*denom);
        numgamma = 16.0*((m(i) - m(k))*var(j) + (m(i) - m(j))*var(k)) +
          4.0*w2*(2.0*m(i) - m(j) - m(k));
        gamma_ijk = -epsilon_ijk*var(i)*numgamma/denom;
        epsilon_nn(i) = epsilon_nn(i) + epsilon_ijk*n(j)*n(k);
        gamma_nn(i) = gamma_nn(i) + gamma_ijk*n(j)*n(k);
      }
    }
  }
  for (i = 0; i < S; i++) {
    dvdt(i)   = n(i)  * (b(i) - alpha_n(i) - epsilon_nn(i)) * smoothstep(n(i)*1.0e6);
    dvdt(i+S) = h2(i) * (g(i) - beta_n(i)  - gamma_nn(i))   * smoothstep(n(i)*1.0e5);
  }
  return(List::create(dvdt));
}


// [[Rcpp::export]]
List hierHOI(double time, NumericVector state, List pars) {
  int S = pars["S"], i, j, k;
  double w = pars["w"], w2 = w*w, W = pars["W"], W2 = W*W, kappa = pars["kappa"];
  double z0 = pars["z0"], sv, svH, dm, ddm, epsilon_ijk, gamma_ijk;
  NumericVector dvdt(2*S), n = state[Range(0, S - 1)], m = state[Range(S, 2*S - 1)];
  NumericVector sig = pars["sigma"], h2 = pars["h2"], var(S);
  NumericMatrix alpha(S,S), beta(S,S);
  List bg = b_hier(sig, pars["theta"], m);
  NumericVector b = bg[0], g = bg[1];
  NumericVector alpha_n(S), beta_n(S), epsilon_nn(S), gamma_nn(S);
  for (i = 0; i < S; i++) var(i) = sig(i)*sig(i);
  for (i = 0; i < S; i++) {
    alpha_n(i) = 0.0;
    beta_n(i) = 0.0;
    epsilon_nn(i) = 0.0;
    gamma_nn(i) = 0.0;
    for (j = 0; j < S; j++) {
      dm = m(i) - m(j);
      sv = 2.0*(var(i) + var(j)) + w2;
      alpha(i,j) = sigmoid(dm / sqrt(sv));
      beta(i,j) = var(i)*exp(-dm*dm / sv) / sqrt(sv*M_PI);
      alpha_n(i) = alpha_n(i) + alpha(i,j)*n(j);
      beta_n(i) = beta_n(i) + beta(i,j)*n(j);
      for (k = 0; k < S; k++) {
        svH = W2 + 2.0*var(i) + var(j)/2.0 + var(k)/2.0;
        ddm = z0 + m(i) - m(j)/2.0 - m(k)/2.0;
        epsilon_ijk = kappa * sigmoid(ddm / sqrt(svH));
        gamma_ijk = kappa * var(i) * exp(-ddm*ddm / svH) / sqrt(svH*M_PI);
        epsilon_nn(i) = epsilon_nn(i) + epsilon_ijk*n(j)*n(k);
        gamma_nn(i) = gamma_nn(i) + gamma_ijk*n(j)*n(k);
      }
    }
  }
  for (i = 0; i < S; i++) {
    dvdt(i)   = n(i)  * (b(i) - alpha_n(i) - epsilon_nn(i)) * smoothstep(n(i)*1.0e6);
    dvdt(i+S) = h2(i) * (g(i) - beta_n(i)  - gamma_nn(i))   * smoothstep(n(i)*1.0e5);
  }
  return(List::create(dvdt));
}
