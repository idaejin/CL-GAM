#include <TMB.hpp>

/* Composite-link Poisson CL-GAMM (Case A style):
 *   eta = X * beta + Z * u
 *   gamma = e * exp(eta)
 *   mu = C * gamma   (partition: accumulate by groups)
 *   y ~ Poisson(mu)
 *   u ~ N(0, G), with diagonal precision
 *     Ginv = G1inv_n / tau2[0] + G2inv_n / tau2[1] + Gnl_n / tau2[2]
 *   matching pois_SOP anisotropic spatial + one nl smooth.
 */
template<class Type>
Type objective_function<Type>::operator()()
{
  using namespace density;

  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_MATRIX(Z);
  DATA_VECTOR(e);
  DATA_IVECTOR(groups);   // 0-based fine -> coarse
  DATA_INTEGER(n_coarse);
  DATA_VECTOR(G1inv_n);   // length = ncol(Z)
  DATA_VECTOR(G2inv_n);
  DATA_VECTOR(Gnl_n);     // 0 if no nl smooth

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(u);
  PARAMETER_VECTOR(log_tau2); // length 2 (spatial) or 3 (+ nl)

  int n_fine = X.rows();
  int q = u.size();
  int n_tau = log_tau2.size();

  Type nll = Type(0);

  // Diagonal precision for random effects
  vector<Type> tau2(n_tau);
  for (int k = 0; k < n_tau; k++) tau2(k) = exp(log_tau2(k));

  for (int j = 0; j < q; j++) {
    Type prec = G1inv_n(j) / tau2(0) + G2inv_n(j) / tau2(1);
    if (n_tau > 2) prec += Gnl_n(j) / tau2(2);
    // Floor for numerical safety
    if (prec < Type(1e-12)) prec = Type(1e-12);
    nll -= dnorm(u(j), Type(0), Type(1) / sqrt(prec), true);
  }

  vector<Type> eta = X * beta + Z * u;
  vector<Type> gamma(n_fine);
  for (int i = 0; i < n_fine; i++) {
    // Soft clip extreme linear predictors
    Type eti = eta(i);
    if (eti > Type(20)) eti = Type(20);
    if (eti < Type(-20)) eti = Type(-20);
    gamma(i) = e(i) * exp(eti);
    eta(i) = eti;
  }

  vector<Type> mu(n_coarse);
  mu.setZero();
  for (int i = 0; i < n_fine; i++) {
    mu(groups(i)) += gamma(i);
  }

  for (int c = 0; c < n_coarse; c++) {
    Type muc = mu(c);
    if (muc < Type(1e-12)) muc = Type(1e-12);
    nll -= dpois(y(c), muc, true);
  }

  ADREPORT(eta);
  ADREPORT(mu);
  REPORT(eta);
  REPORT(mu);
  REPORT(gamma);
  return nll;
}
