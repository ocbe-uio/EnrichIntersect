
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <limits>


// Convert R 1-based integer indices to Armadillo 0-based index vector.
// Uses int / unsigned int explicitly, no arma::uword in user code.
arma::uvec r_index_to_arma_uvec_original(
    const Rcpp::IntegerVector& idx,
    const unsigned int n
) {
  std::vector<unsigned int> tmp;
  tmp.reserve(static_cast<unsigned int>(idx.size()));

  for (int j = 0; j < idx.size(); ++j) {
    if (idx[j] == NA_INTEGER || idx[j] < 1) {
      continue;
    }

    const unsigned int idx0 = static_cast<unsigned int>(idx[j] - 1);

    if (idx0 < n) {
      tmp.push_back(idx0);
    }
  }

  arma::uvec out(tmp.size());

  for (unsigned int i = 0; i < tmp.size(); ++i) {
    out[i] = tmp[i];
  }

  return out;
}


// Original one-sided KS-like enrichment score:
//
//   ES = max(F1 - F2)
//
// ordered_scores:
//   scores sorted decreasingly, as in the original R code:
//
//     myList <- sort(myList, decreasing = TRUE)
//
// hit_mask:
//   0/1 vector in the same ranked order as ordered_scores.
double calc_es_original(
    const arma::vec& ordered_scores,
    const arma::uvec& hit_mask,
    const double alpha
) {
  const unsigned int n = static_cast<unsigned int>(ordered_scores.n_elem);
  const unsigned int k = static_cast<unsigned int>(arma::accu(hit_mask));

  if (n == 0 || k == 0 || k >= n) {
    return NA_REAL;
  }

  arma::vec hit_numeric = arma::conv_to<arma::vec>::from(hit_mask);
  arma::vec miss_numeric = 1.0 - hit_numeric;

  arma::vec weights =
    arma::pow(arma::abs(ordered_scores), alpha) % hit_numeric;

  const double weight_sum = arma::accu(weights);

  if (weight_sum <= 0.0 || !std::isfinite(weight_sum)) {
    return NA_REAL;
  }

  arma::vec F1 = arma::cumsum(weights) / weight_sum;
  arma::vec F2 =
    arma::cumsum(miss_numeric) / static_cast<double>(n - k);

  arma::vec diff = F1 - F2;

  return diff.max();
}


// Partial Fisher-Yates sample without replacement.
// Returns k unique positions from 0:(n - 1).
//
// This avoids shuffling all n positions when only k are needed.
arma::uvec sample_positions_partial(
    const unsigned int n,
    const unsigned int k
) {
  arma::uvec pool(n);

  for (unsigned int i = 0; i < n; ++i) {
    pool[i] = i;
  }

  for (unsigned int i = 0; i < k; ++i) {
    const unsigned int j = i + static_cast<unsigned int>(
      std::floor(R::runif(0.0, static_cast<double>(n - i)))
    );

    const unsigned int tmp = pool[i];
    pool[i] = pool[j];
    pool[j] = tmp;
  }

  return pool.head(k);
}


// [[Rcpp::export]]
Rcpp::List enrichment_core_original(
    const arma::mat& x,
    const Rcpp::List& set_indices,
    const double alpha = 0.0,
    const bool normalize = true,
    const int permute_n = 100
) {
  if (permute_n < 0) {
    Rcpp::stop("'permute_n' must be non-negative.");
  }

  // I decided to use simple 'unsigned int' rather than 'arma::uword' but sacrifice length
  if (x.n_rows > std::numeric_limits<unsigned int>::max()) {
    Rcpp::stop("Too many rows in 'x' for unsigned int indexing.");
  }

  if (x.n_cols > std::numeric_limits<unsigned int>::max()) {
    Rcpp::stop("Too many columns in 'x' for unsigned int indexing.");
  }

  if (set_indices.size() > std::numeric_limits<unsigned int>::max()) {
    Rcpp::stop("Too many custom sets for unsigned int indexing.");
  }

  const unsigned int n_features =
    static_cast<unsigned int>(x.n_rows);

  const unsigned int n_cols =
    static_cast<unsigned int>(x.n_cols);

  const unsigned int n_groups =
    static_cast<unsigned int>(set_indices.size());

  const unsigned int permute_n_u =
    static_cast<unsigned int>(permute_n);

  arma::mat S(
    n_cols,
    n_groups,
    arma::fill::value(NA_REAL)
  );

  arma::mat pvalue(
    n_cols,
    n_groups,
    arma::fill::value(NA_REAL)
  );

  // Convert all R 1-based set indices to Armadillo 0-based index vectors.
  std::vector<arma::uvec> sets;
  sets.reserve(n_groups);

  for (unsigned int g = 0; g < n_groups; ++g) {
    Rcpp::IntegerVector idx_r = set_indices[g];

    sets.push_back(
      r_index_to_arma_uvec_original(
        idx_r,
        n_features
      )
    );
  }

  Rcpp::RNGScope scope;

  for (unsigned int col = 0; col < n_cols; ++col) {
    arma::vec scores = x.col(col);

    // R's sort() removes NA by default.
    // Here we remove non-finite values.
    arma::uvec finite_idx = arma::find_finite(scores);

    if (finite_idx.n_elem == 0) {
      continue;
    }

    arma::vec finite_scores = scores.elem(finite_idx);

    // Original R code:
    //
    //   myList <- sort(myList, decreasing = TRUE)
    //
    arma::uvec order_local =
      arma::sort_index(finite_scores, "descend");

    arma::uvec ordered_idx =
      finite_idx.elem(order_local);

    arma::vec ordered_scores =
      scores.elem(ordered_idx);

    const unsigned int n_ranked =
      static_cast<unsigned int>(ordered_idx.n_elem);

    for (unsigned int g = 0; g < n_groups; ++g) {
      const arma::uvec& set_idx = sets[g];

      if (set_idx.n_elem == 0) {
        S(col, g) = NA_REAL;
        pvalue(col, g) = NA_REAL;
        continue;
      }

      arma::uvec global_hit_mask(
        n_features,
        arma::fill::zeros
      );

      global_hit_mask.elem(set_idx).ones();

      arma::uvec hit_mask =
        global_hit_mask.elem(ordered_idx);

      const unsigned int k =
        static_cast<unsigned int>(arma::accu(hit_mask));

      if (k == 0 || k >= n_ranked) {
        S(col, g) = NA_REAL;
        pvalue(col, g) = NA_REAL;
        continue;
      }

      const double observed = calc_es_original(
        ordered_scores,
        hit_mask,
        alpha
      );

      if (!std::isfinite(observed)) {
        S(col, g) = NA_REAL;
        pvalue(col, g) = NA_REAL;
        continue;
      }

      arma::vec permute_S(
        permute_n_u,
        arma::fill::value(NA_REAL)
      );

      for (unsigned int b = 0; b < permute_n_u; ++b) {
        arma::uvec perm_hit_mask(
          n_ranked,
          arma::fill::zeros
        );

        arma::uvec sampled_pos =
          sample_positions_partial(
            n_ranked,
            k
          );

        perm_hit_mask.elem(sampled_pos).ones();

        permute_S[b] = calc_es_original(
          ordered_scores,
          perm_hit_mask,
          alpha
        );
      }

      arma::uvec valid_perm =
        arma::find_finite(permute_S);

      if (valid_perm.n_elem == 0) {
        S(col, g) = observed;
        pvalue(col, g) = NA_REAL;
        continue;
      }

      double observed_for_p = observed;
      arma::vec permute_for_p = permute_S;

      if (normalize) {
        arma::uvec pos_perm =
          arma::find(permute_S >= 0.0);

        if (pos_perm.n_elem > 0) {
          const double denom =
            arma::mean(permute_S.elem(pos_perm));

          if (std::isfinite(denom) && denom != 0.0) {
            observed_for_p = observed / denom;
            permute_for_p = permute_S / denom;
          }
        }
      }

      S(col, g) = observed_for_p;

      arma::uvec valid_norm =
        arma::find_finite(permute_for_p);

      if (valid_norm.n_elem == 0) {
        pvalue(col, g) = NA_REAL;
      } else {
        arma::vec valid_values =
          permute_for_p.elem(valid_norm);

        arma::uvec extreme =
          arma::find(valid_values >= observed_for_p);

        // Original R code used no +1 correction:
        //
        //   sum(permute_S >= S) / valid permutations
        //
        pvalue(col, g) =
          static_cast<double>(extreme.n_elem) /
          static_cast<double>(valid_norm.n_elem);
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("S") = S,
    Rcpp::Named("pvalue") = pvalue
  );
}