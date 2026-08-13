#include <TMB.hpp>

/* Composite-link Poisson CL-GAMM (Laplace):
 *   eta = X * beta + Z * u
 *   gamma = e * exp(eta)
 *   mu_c accumulated by partition groups
 *   y ~ Poisson(mu)
 *   u ~ N(0, G) with diagonal precision
 *     Ginv(j) = sum_k Ginv_n(j, k) / tau2(k)
 *   Columns of Ginv_n are spatial tau1, tau2, then each smooth, then
 *   optional coarse iid RE (last column).
 */
template<class Type>
Type objective_function<Type>::operator()()
{
  using namespace density;

  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_MATRIX(Z);
  DATA_VECTOR(e);
  DATA_IVECTOR(groups);
  DATA_INTEGER(n_coarse);
  DATA_MATRIX(Ginv_n);
  DATA_INTEGER(n_re);

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(u);
  PARAMETER_VECTOR(log_tau2);

  int n_fine = X.rows();
  int q = u.size();
  int n_tau = log_tau2.size();

  Type nll = Type(0);

  vector<Type> tau2(n_tau);
  for (int k = 0; k < n_tau; k++) tau2(k) = exp(log_tau2(k));

  for (int j = 0; j < q; j++) {
    Type prec = Type(0);
    for (int k = 0; k < n_tau; k++) {
      prec += Ginv_n(j, k) / tau2(k);
    }
    if (prec < Type(1e-12)) prec = Type(1e-12);
    nll -= dnorm(u(j), Type(0), Type(1) / sqrt(prec), true);
  }

  vector<Type> eta = X * beta + Z * u;
  for (int i = 0; i < n_fine; i++) {
    if (eta(i) > Type(20)) eta(i) = Type(20);
    if (eta(i) < Type(-20)) eta(i) = Type(-20);
  }

  vector<Type> eta_struct = eta;
  if (n_re > 0) {
    int q_struct = q - n_re;
    for (int i = 0; i < n_fine; i++) {
      Type re_i = Type(0);
      for (int r = 0; r < n_re; r++) {
        re_i += Z(i, q_struct + r) * u(q_struct + r);
      }
      eta_struct(i) = eta(i) - re_i;
    }
  }

  vector<Type> gamma(n_fine);
  for (int i = 0; i < n_fine; i++) {
    gamma(i) = e(i) * exp(eta(i));
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

  ADREPORT(eta_struct);
  ADREPORT(eta);
  ADREPORT(mu);
  REPORT(eta_struct);
  REPORT(eta);
  REPORT(mu);
  REPORT(gamma);
  return nll;
}
