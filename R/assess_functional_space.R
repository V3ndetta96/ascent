#' Assess Functional Space Quality
#'
#' @description
#' Diagnostic tool to evaluate the quality of the PCoA-based functional space
#' constructed from a trait matrix. Intended for use before running core
#' analyses (\code{asc_paired}, \code{asc_pairwise}, etc.) to identify
#' potential distortion from negative eigenvalues or excessive dimensionality
#' reduction.
#'
#' @param traits A data frame or matrix of functional traits.
#' @param dist_method Distance metric for the trait matrix. Default is \code{"gower"}.
#' @param dim_retention Character. Method for dimensionality reduction.
#'   Options are \code{"variance"} or \code{"broken_stick"}. Default is \code{"variance"}.
#' @param var_tol Numeric. Proportion of cumulative variance to retain (used
#'   when \code{dim_retention = "variance"}). Default is 0.80.
#'
#' @return A list of class \code{ascfcd_space} with:
#' \describe{
#'   \item{quality}{Numeric. PCoA quality index (\code{sum(positive_eig) / sum(abs(all_eig))}). Values above 0.80 are generally acceptable.}
#'   \item{n_species}{Integer. Number of species in the trait matrix.}
#'   \item{n_traits}{Integer. Number of traits.}
#'   \item{k_retained}{Integer. Number of PCoA axes retained.}
#'   \item{var_retained}{Numeric. Cumulative variance explained by retained axes.}
#'   \item{axis_var}{Numeric vector. Relative variance per retained axis.}
#'   \item{n_neg_eigenvalues}{Integer. Count of negative eigenvalues.}
#'   \item{neg_eigenvalue_pct}{Numeric. Percentage of total absolute eigenvalue represented by negative eigenvalues.}
#'   \item{eigenvalues}{Numeric vector. All raw eigenvalues from PCoA.}
#' }
#'
#' @details
#' The Gower distance is not strictly Euclidean, which means PCoA may produce
#' negative eigenvalues. This function reports the full eigenvalue spectrum
#' and computes a quality index so users can decide whether the functional
#' space is an adequate representation of the original distance matrix.
#'
#' A quality below 0.80 indicates that negative eigenvalues represent more
#' than 20\% of the total variation. In such cases, consider:
#' \itemize{
#'   \item Reviewing the trait data for problematic variables.
#'   \item Applying a Cailliez or Lingoes correction to the distance matrix before ordination.
#'   \item Using a different distance metric (e.g., \code{"euclidean"} for purely continuous traits).
#' }
#'
#' @examples
#' traits <- data.frame(
#'   Mass = c(15, 30, 60, 150, 400),
#'   Beak = c(10, 15, 28,  45,  85),
#'   Diet = factor(c(0, 1, 1, 1, 0))
#' )
#' rownames(traits) <- paste0("Sp", 1:5)
#'
#' diag <- assess_functional_space(traits)
#' diag
#'
#' @export
assess_functional_space <- function(traits, dist_method = "gower",
                                     dim_retention = c("variance", "broken_stick"),
                                     var_tol = 0.80) {

  dim_retention <- match.arg(dim_retention)

  # Compute distance and PCoA
  d_mat <- cluster::daisy(traits, metric = dist_method)
  pcoa_res <- stats::cmdscale(d_mat, k = nrow(traits) - 1, eig = TRUE)

  eig_raw <- pcoa_res$eig
  eig_pos <- eig_raw[eig_raw > 1e-8]
  n_neg <- sum(eig_raw < -1e-8)
  neg_pct <- if (sum(abs(eig_raw)) > 0) {
    sum(abs(eig_raw[eig_raw < -1e-8])) / sum(abs(eig_raw)) * 100
  } else {
    0
  }

  quality <- sum(eig_pos) / sum(abs(eig_raw))

  # Axis retention (same logic as .build_functional_space)
  rel_var <- eig_pos / sum(eig_pos)
  cum_var <- cumsum(rel_var)

  if (dim_retention == "variance") {
    k_valid <- min(which(cum_var >= var_tol))
  } else if (dim_retention == "broken_stick") {
    n_eig <- length(eig_pos)
    bs_expected <- sapply(1:n_eig, function(k) sum(1 / (k:n_eig)) / n_eig)
    k_valid <- sum(rel_var > bs_expected)
  }
  if (k_valid < 2) k_valid <- 2

  output <- list(
    quality            = quality,
    n_species          = nrow(traits),
    n_traits           = ncol(traits),
    k_retained         = k_valid,
    var_retained       = cum_var[k_valid],
    axis_var           = rel_var[1:k_valid],
    n_neg_eigenvalues  = n_neg,
    neg_eigenvalue_pct = neg_pct,
    eigenvalues        = eig_raw
  )
  class(output) <- "ascfcd_space"
  return(output)
}


#' @export
print.ascfcd_space <- function(x, ...) {
  cat("==================================================\n")
  cat(" Functional Space Diagnostic (ascent)\n")
  cat("==================================================\n\n")
  cat(sprintf("  Species:           %d\n", x$n_species))
  cat(sprintf("  Traits:            %d\n", x$n_traits))
  cat(sprintf("  PCoA Quality:      %.1f%%\n", x$quality * 100))
  cat(sprintf("  Axes Retained:     %d (%.1f%% variance)\n",
              x$k_retained, x$var_retained * 100))
  cat(sprintf("  Negative Eigenvalues: %d (%.1f%% of total)\n",
              x$n_neg_eigenvalues, x$neg_eigenvalue_pct))

  if (x$quality < 0.80) {
    cat("\n  [!] WARNING: Quality < 80%. Consider Cailliez correction\n")
    cat("      or reviewing trait data.\n")
  } else {
    cat("\n  Quality is adequate for analysis.\n")
  }
  invisible(x)
}


#' @export
summary.ascfcd_space <- function(object, ...) {
  print.ascfcd_space(object, ...)
  cat("\n--- Per-Axis Variance ---\n")
  axis_df <- data.frame(
    Axis     = paste0("Axis_", seq_along(object$axis_var)),
    Variance = round(object$axis_var * 100, 2),
    Cumul    = round(cumsum(object$axis_var) / sum(object$axis_var) * object$var_retained * 100, 2)
  )
  print(axis_df, row.names = FALSE)
  cat("\n--- Full Eigenvalue Spectrum ---\n")
  eig <- object$eigenvalues
  cat(sprintf("  Positive: %d | Negative: %d | Total: %d\n",
              sum(eig > 1e-8), sum(eig < -1e-8), length(eig)))
  invisible(object)
}
