// Hot-path kernels for SOP / CLMM (RcppArmadillo)
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Schur SOP solve + ED diagonal (Ayma nonsymmetric system)
//'
//' @param XtX p x p
//' @param ZtX q x p
//' @param ZtZ q x q
//' @param ZtXtZ (p+q) x q
//' @param u length p+q
//' @param G length q (1/Ginv)
//' @param A11inv_cached optional p x p inverse of XtX (empty = compute)
//' @return list b_fixed, b_random, dZtNZ, A11inv
//' @keywords internal
// [[Rcpp::export]]
List sop_solve_schur_cpp(const arma::mat& XtX,
                         const arma::mat& ZtX,
                         const arma::mat& ZtZ,
                         const arma::mat& ZtXtZ,
                         const arma::vec& u,
                         const arma::vec& G,
                         const arma::mat& A11inv_cached) {
  const uword p = XtX.n_cols;
  const uword q = G.n_elem;

  arma::mat A11inv;
  if (A11inv_cached.n_elem == 0) {
    if (!inv_sympd(A11inv, XtX)) {
      A11inv = pinv(XtX);
    }
  } else {
    A11inv = A11inv_cached;
  }

  // A12 = t(diag(G) * ZtX); A21 = ZtX; A22 = t(diag(G) * ZtZ) + I
  arma::mat GZ = ZtX.each_col() % G;   // diag(G) %*% ZtX
  arma::mat A12 = GZ.t();              // p x q
  const arma::mat& A21 = ZtX;          // q x p
  arma::mat A22 = (ZtZ.each_col() % G).t();
  A22.diag() += 1.0;

  arma::vec u1 = u.head(p);
  arma::vec u2 = u.tail(q);

  arma::mat A11inv_A12 = A11inv * A12;
  arma::mat S = A22 - A21 * A11inv_A12;
  arma::vec rhs2 = u2 - A21 * (A11inv * u1);

  arma::mat Sinv;
  arma::vec b2;
  arma::mat I = eye(q, q);
  if (!solve(Sinv, S, I, solve_opts::likely_sympd + solve_opts::no_approx)) {
    Sinv = pinv(S);
  }
  b2 = Sinv * rhs2;

  arma::vec b1 = A11inv * (u1 - A12 * b2);

  // H22 = Sinv, H21 = -Sinv * A21 * A11inv
  arma::mat H21 = -Sinv * (A21 * A11inv);
  arma::mat H_bottom = join_horiz(H21, Sinv); // q x (p+q)

  // dZtNZ[j] = sum_k H_bottom(j,k) * ZtXtZ(k,j)
  arma::vec dZtNZ(q);
  for (uword j = 0; j < q; ++j) {
    dZtNZ(j) = dot(H_bottom.row(j), ZtXtZ.col(j));
  }

  arma::vec b_random = G % b2;

  return List::create(
    Named("b.fixed") = b1,
    Named("b.random") = b_random,
    Named("dZtNZ") = dZtNZ,
    Named("A11inv") = A11inv
  );
}

//' Partition (rowsum) aggregation: C %*% (gamma * A) for 0-1 C
//' @param gamma length n_fine
//' @param A n_fine x p
//' @param groups 1-based group id length n_fine
//' @param n_coarse number of coarse cells
//' @keywords internal
// [[Rcpp::export]]
arma::mat comp_mul_groups_cpp(const arma::vec& gamma,
                              const arma::mat& A,
                              const arma::uvec& groups,
                              const int n_coarse) {
  const uword n = A.n_rows;
  const uword p = A.n_cols;
  arma::mat out(n_coarse, p, fill::zeros);
  for (uword i = 0; i < n; ++i) {
    uword g = groups(i) - 1; // 1-based -> 0-based
    out.row(g) += gamma(i) * A.row(i);
  }
  return out;
}

//' Working cross-products for CLMM given CGX, CGZ (coarse x p / q)
//' @keywords internal
// [[Rcpp::export]]
List clmm_crossprod_cpp(const arma::mat& CGX,
                        const arma::mat& CGZ,
                        const arma::vec& z,
                        const arma::vec& w) {
  arma::vec sw = sqrt(1.0 / w);
  arma::mat XtX = (CGX.each_col() % sw).t() * (CGX.each_col() % sw);
  arma::mat ZtZ = (CGZ.each_col() % sw).t() * (CGZ.each_col() % sw);
  arma::mat XtZ = CGX.t() * (CGZ.each_col() / w);
  arma::mat ZtX = XtZ.t();
  arma::vec Xty = CGX.t() * z;
  arma::vec Zty = CGZ.t() * z;
  double yty = accu((z % z) % w);
  arma::mat ZtXtZ = join_vert(XtZ, ZtZ);
  arma::vec u = join_vert(Xty, Zty);
  return List::create(
    Named("XtX") = XtX,
    Named("XtZ") = XtZ,
    Named("ZtX") = ZtX,
    Named("ZtZ") = ZtZ,
    Named("Xty") = Xty,
    Named("Zty") = Zty,
    Named("yty") = yty,
    Named("ZtXtZ") = ZtXtZ,
    Named("u") = u
  );
}

//' Form B' diag(w) B + apply (dense) for Camarda update RHS B'(w*z)
//' @keywords internal
// [[Rcpp::export]]
List btWb_cpp(const arma::mat& B,
              const arma::vec& w,
              const arma::vec& z) {
  arma::vec sw = sqrt(w);
  arma::mat BtWB = (B.each_col() % sw).t() * (B.each_col() % sw);
  arma::vec rhs = B.t() * (w % z);
  return List::create(Named("BtWB") = BtWB, Named("rhs") = rhs);
}
