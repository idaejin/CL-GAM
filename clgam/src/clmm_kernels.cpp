// Hot-path kernels for SOP / CLMM (RcppArmadillo)
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Schur SOP solve + ED diagonal (Ayma nonsymmetric system)
//'
//' Algebraic note (verified numerically, see clgam SOP review notes): with
//' A12 = XtZ \%*\% diag(G), A22 = ZtZ \%*\% diag(G) + I (both G-dependent), and
//' A11inv fixed,
//'   S = A22 - ZtX \%*\% A11inv \%*\% A12 = N \%*\% diag(G) + I,
//'   rhs2 = u2 - ZtX \%*\% A11inv \%*\% u1,
//' where N = ZtZ - ZtX \%*\% A11inv \%*\% XtZ and rhs2 do NOT depend on G at
//' all. Within one outer PIRLS iteration only G changes across inner SOP
//' iterations (XtX, ZtX, ZtZ, u come from the working system formed once
//' before the inner loop), so N and rhs2 are computed once (when the cache
//' is empty) and reused: forming S per inner iteration drops from two
//' O(q^2 p) matrix products to one O(q^2) column rescale.
//'
//' Also note S = N \%*\% diag(G) + I is NOT symmetric in general (only if G
//' is constant), even though N itself is symmetric -- confirmed
//' numerically (||S - t(S)|| was large, not zero, on random test inputs).
//' The previous version solved it via `solve_opts::likely_sympd`, which
//' assumes symmetry; this version uses a general (LU-based) inverse, which
//' is correct regardless of symmetry and matches what the pure-R fallback
//' (`.sop_solve_schur_R`, plain `solve()`) has always done.
//'
//' @param XtX p x p
//' @param ZtX q x p
//' @param ZtZ q x q
//' @param ZtXtZ (p+q) x q
//' @param u length p+q
//' @param G length q (1/Ginv)
//' @param A11inv_cached optional p x p inverse of XtX (empty = compute)
//' @param N_cached optional q x q G-free Schur complement (empty = compute)
//' @param rhs2_cached optional length-q G-free RHS (empty = compute)
//' @return list b_fixed, b_random, dZtNZ, A11inv, N, rhs2
//' @keywords internal
// [[Rcpp::export]]
List sop_solve_schur_cpp(const arma::mat& XtX,
                         const arma::mat& ZtX,
                         const arma::mat& ZtZ,
                         const arma::mat& ZtXtZ,
                         const arma::vec& u,
                         const arma::vec& G,
                         const arma::mat& A11inv_cached,
                         const arma::mat& N_cached,
                         const arma::vec& rhs2_cached) {
  const uword p = XtX.n_cols;
  const uword q = G.n_elem;

  arma::vec u1 = u.head(p);
  arma::vec u2 = u.tail(q);

  const bool have_cache = (A11inv_cached.n_elem > 0) &&
                           (N_cached.n_elem > 0) &&
                           (rhs2_cached.n_elem > 0);

  arma::mat A11inv;
  arma::mat N;
  arma::vec rhs2;
  arma::mat XtZ = ZtX.t();

  if (!have_cache) {
    if (!inv_sympd(A11inv, XtX)) {
      A11inv = pinv(XtX);
    }
    N = ZtZ - ZtX * (A11inv * XtZ);
    rhs2 = u2 - ZtX * (A11inv * u1);
  } else {
    A11inv = A11inv_cached;
    N = N_cached;
    rhs2 = rhs2_cached;
  }

  // S = N %*% diag(G) + I (column-wise rescale of the cached, G-free N).
  arma::mat S = N.each_row() % G.t();
  S.diag() += 1.0;

  // A12 = XtZ %*% diag(G) (p x q), still needed for b1; cheap (O(p*q)).
  arma::mat A12 = XtZ.each_row() % G.t();

  arma::mat Sinv;
  if (!inv(Sinv, S)) {
    Sinv = pinv(S);
  }
  arma::vec b2 = Sinv * rhs2;

  arma::vec b1 = A11inv * (u1 - A12 * b2);

  // H22 = Sinv, H21 = -Sinv * ZtX * A11inv
  arma::mat H21 = -Sinv * (ZtX * A11inv);
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
    Named("A11inv") = A11inv,
    Named("N") = N,
    Named("rhs2") = rhs2
  );
}

//' Partition (rowsum) aggregation: C \%*\% (gamma * A) for 0-1 C
//'
//' Compiled helper for 0-1 composition matrices.
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
