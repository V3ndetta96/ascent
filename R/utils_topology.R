# ==============================================================================
# INTERNAL UTILITY FUNCTIONS — SHARED ENGINES
# ==============================================================================
# These are hidden (.dot) functions that serve as the single source of truth
# for the core mathematical engines used across the package.

#' Build the Functional PCoA Space from a Trait Matrix
#'
#' @description
#' Internal function. Computes Gower (or other) distance, runs PCoA,
#' selects axes via variance threshold or broken stick, and reports quality.
#'
#' @param traits Data frame or matrix of functional traits.
#' @param dist_method Distance metric (default \code{"gower"}).
#' @param dim_retention Method for axis selection (\code{"variance"} or \code{"broken_stick"}).
#' @param var_tol Proportion of cumulative variance to retain.
#'
#' @return A list with: \code{F_mat}, \code{k_retained}, \code{var_retained},
#'   \code{axis_var}, \code{quality}.
#'
#' @details
#' **Negative eigenvalues:** Gower distance is not strictly Euclidean, so
#' \code{cmdscale} may produce negative eigenvalues. This function filters
#' them out and reports a quality metric (\code{quality = sum(positive_eig) /
#' sum(abs(all_eig))}). If quality falls below 0.80, a warning is issued.
#' No Cailliez or Lingoes correction is applied; the distances in the
#' retained PCoA space are an approximation of the original Gower distances.
#' For matrices with high proportions of binary/categorical traits, the
#' distortion may be substantial.
#'
#' @noRd
.build_functional_space <- function(traits, dist_method = "gower",
                                     dim_retention = "variance", var_tol = 0.80) {

  d_mat <- cluster::daisy(traits, metric = dist_method)
  pcoa_res <- stats::cmdscale(d_mat, k = nrow(traits) - 1, eig = TRUE)

  eig_raw <- pcoa_res$eig
  eig_pos <- eig_raw[eig_raw > 1e-8]
  rel_var <- eig_pos / sum(eig_pos)
  cum_var <- cumsum(rel_var)

  # PCoA quality: proportion of total absolute eigenvalue explained by positive axes

  quality <- sum(eig_pos) / sum(abs(eig_raw))
  if (quality < 0.80) {
    warning(sprintf(
      paste0("PCoA quality is low (%.1f%%). Negative eigenvalues represent ",
             "%.1f%% of total variation. Consider applying Cailliez correction ",
             "or reviewing trait data."),
      quality * 100, (1 - quality) * 100
    ))
  }

  if (dim_retention == "variance") {
    k_valid <- min(which(cum_var >= var_tol))
  } else if (dim_retention == "broken_stick") {
    n_eig <- length(eig_pos)
    bs_expected <- sapply(1:n_eig, function(k) sum(1 / (k:n_eig)) / n_eig)
    k_valid <- sum(rel_var > bs_expected)
  }
  if (k_valid < 2) k_valid <- 2

  axis_var <- rel_var[1:k_valid]
  total_var_retained <- cum_var[k_valid]
  F_mat <- pcoa_res$points[, 1:k_valid, drop = FALSE]
  colnames(F_mat) <- paste0("Axis_", 1:k_valid)

  list(
    F_mat = F_mat,
    k_retained = k_valid,
    var_retained = total_var_retained,
    axis_var = axis_var,
    quality = quality
  )
}


#' Compute Core Topology (FDis, FRic) for a Single Community
#'
#' @description
#' Internal function. Calculates abundance-weighted functional dispersion
#' (Laliberte & Legendre 2010) and convex hull volume for a single
#' community vector in the retained functional space.
#'
#' @param abund_vec Numeric vector of (relative) abundances.
#' @param F_space Matrix of species coordinates in PCoA space.
#'
#' @return A list with \code{FDis} and \code{FRic}.
#'   \code{FRic} is \code{NA_real_} when the convex hull is not computable.
#' @noRd
.calc_core_topology <- function(abund_vec, F_space) {
  idx <- abund_vec > 0
  n_pres <- sum(idx)

  if (n_pres < 2) return(list(FDis = 0, FRic = NA_real_))

  p_pres <- abund_vec[idx] / sum(abund_vec[idx])
  F_pres <- F_space[idx, , drop = FALSE]

  centroid <- colSums(p_pres * F_pres)
  dist_to_c <- sqrt(rowSums(sweep(F_pres, 2, centroid, "-")^2))
  fdis <- sum(p_pres * dist_to_c)

  fric <- NA_real_
  if (n_pres > ncol(F_space)) {
    fric <- tryCatch(
      geometry::convhulln(F_pres, options = "FA")$vol,
      error = function(e) NA_real_
    )
  }

  list(FDis = fdis, FRic = fric)
}
