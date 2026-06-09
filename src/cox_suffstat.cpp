// cox_suffstat.cpp
// Single-pass Breslow accumulation for Cox fine-mapping sufficient statistics.
// Produces the pieces needed to build XtX = A - B'B and Xty = X' M,
// without ever forming an n x n matrix.
//
// Breslow ties handled by two-stage block flattening: all individuals sharing
// the same event/censoring time share the risk-set sums S0 / S1 evaluated at
// the end of their tied block.

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

// Returns:
//   a      (n)   per-individual cumulative weight  Lambda0(t_i) * exp(eta_i)
//   B      (d x p) risk-set weighted means xbar_j at each EVENT time (Breslow)
//   Xty    (p)   X' M, the Cox score vector
//   d      int   number of events
//
// Inputs assumed: time, status, eta are length n; X is n x p.
// Sorting is done internally (descending time). Ties in time are flattened.
// n_threads: the row scan is sequential; this only sets the OpenMP thread
//   count seen by downstream BLAS-backed crossprods. Default 1.
// [[Rcpp::export]]
Rcpp::List cox_suffstat(const arma::mat& X,
                        const arma::vec& eta,
                        const arma::vec& time,
                        const arma::ivec& status,
                        int n_threads = 1) {

  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;

#ifdef _OPENMP
  if (n_threads < 1) n_threads = 1;
  omp_set_num_threads(n_threads);
#endif

  // descending sort by time; risk set then grows monotonically as we advance
  arma::uvec ord = arma::sort_index(time, "descend");

  // eta shift: subtract max(eta) before exp for numerical stability.
  // The shift cancels in every risk-set ratio (S1/S0, dev/S0) and in a_i
  // (exp(eta_i) * Lambda0), so it leaves all outputs unchanged. Using max(eta)
  // makes the largest weight exactly 1 (no overflow, no over-shrinking).
  double eta_max = eta.max();
  arma::vec  w  = arma::exp(eta - eta_max);  // shifted weights, original order
  arma::vec  a(n, arma::fill::zeros);        // output cumulative weight, original order

  // running risk-set sums
  double         S0 = 0.0;
  arma::rowvec   S1(p, arma::fill::zeros);

  // event bookkeeping, recorded per EVENT (collected latest-first in descending pass)
  std::vector<arma::rowvec> ev_xbar;       // flattened xbar_j (length p), one per event
  // event-time BLOCK bookkeeping for Lambda0 (one entry per event time, not per event)
  std::vector<double> blk_time;            // event-block time
  std::vector<double> blk_inc;             // Breslow increment dev / S0_flat

  arma::uword i = 0;
  while (i < n) {
    // identify tied block [i, j) with identical time
    arma::uword j = i;
    double t_block = time(ord(i));
    while (j < n && time(ord(j)) == t_block) ++j;

    // ----- stage 1: accumulate the whole block into S0, S1 -----
    // No inner OpenMP here: the row scan is inherently sequential (cumsum
    // dependence), and a per-row fork-join would fire ~n times on continuous
    // time, dwarfing any benefit. Armadillo's lazy row add triggers BLAS/SIMD
    // with zero temporary copy. Parallelism is reserved for the crossprods
    // (A, B'B) on the R side via multithreaded BLAS.
    for (arma::uword r = i; r < j; ++r) {
      arma::uword idx = ord(r);
      double wr = w(idx);
      S0 += wr;
      S1 += wr * X.row(idx);
    }

    // ----- stage 2: block-flattened risk-set mean shared by all tied rows -----
    arma::rowvec xbar = S1 / S0;           // Breslow: same for every event in block

    // count events in block; record one xbar per event, one increment per block
    arma::uword dev = 0;
    for (arma::uword r = i; r < j; ++r)
      if (status(ord(r)) == 1) ++dev;

    if (dev > 0) {
      for (arma::uword e = 0; e < dev; ++e) ev_xbar.push_back(xbar);
      blk_time.push_back(t_block);
      blk_inc.push_back((double)dev / S0);  // Breslow Lambda0 increment for this time
    }

    i = j;
  }

  const arma::uword d = ev_xbar.size();

  // ----- build a_i = exp(eta_i) * Lambda0(t_i) -----
  // Lambda0(t_i) = sum over event-time blocks with t_block <= t_i of (dev / S0_flat).
  // Increments are accumulated PER EVENT-TIME BLOCK (not per row): tied events at the
  // same time contribute a single dev/S0 term, otherwise the increment is double-counted.
  // Single ascending sweep over individuals; blocks were created in descending order,
  // so iterate them in reverse to get ascending event times.
  arma::uvec ord_asc = arma::sort_index(time, "ascend");
  double Lam = 0.0;
  long bptr = (long)blk_time.size() - 1;     // smallest event time is the last block pushed
  for (arma::uword r = 0; r < n; ++r) {
    arma::uword idx = ord_asc(r);
    double ti = time(idx);
    while (bptr >= 0 && blk_time[bptr] <= ti) {
      Lam += blk_inc[bptr];
      --bptr;
    }
    a(idx) = w(idx) * Lam;
  }

  // ----- Xty = X' M, with M_i = status_i - a_i -----
  arma::vec M = arma::conv_to<arma::vec>::from(status) - a;
  arma::vec Xty = X.t() * M;

  // ----- assemble B (d x p) from recorded event xbars -----
  arma::mat B(d, p, arma::fill::none);
  for (arma::uword e = 0; e < d; ++e) {
    B.row(e) = ev_xbar[e];
  }

  return Rcpp::List::create(
    Rcpp::Named("a")   = a,        // n: per-individual cumulative weight
    Rcpp::Named("B")   = B,        // d x p: event-time risk-set means
    Rcpp::Named("Xty") = Xty,      // p: Cox score vector
    Rcpp::Named("d")   = (int)d,   // number of events
    Rcpp::Named("M")   = M         // n: martingale residuals (debug / reuse)
  );
}
