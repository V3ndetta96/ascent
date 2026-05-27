#' Compute Null Models for Functional Centroid Displacement
#'
#' @description
#' Evaluates the statistical significance and magnitude of functional displacement
#' against two null expectations: Neutral Drift (Model 1: Trait shuffling) and
#' Structural Resilience (Model 2: Abundance shuffling).
#'
#' @param x An object of class `ascfcd` (paired) or `ascfcd_pw` (pairwise).
#' @param n_perm Integer. Number of permutations (default: 999).
#' @param seed Integer. Seed for random number generation to ensure reproducibility.
#'
#' @return The original object appended with a `$null_models` data frame containing
#' p-values and Standardized Effect Sizes (SES) for both models.
#' @export
#'
#' @importFrom stats sd
#' @examples
#' \dontrun{
#' res <- asc_paired(traits, ref, comp)
#' res_null <- asc_null(res, n_perm = 999)
#' }
asc_null <- function(x, n_perm = 999, seed = 42) {
  UseMethod("asc_null")
}

#' @export
asc_null.ascfcd <- function(x, n_perm = 999, seed = 42) {
  if (!is.null(seed)) set.seed(seed)

  cat(sprintf("\nRunning %d permutations for paired analysis...\n", n_perm))

  F_mat <- x$F_matrix
  rel_ref <- x$rel_ref
  rel_comp <- x$rel_comp
  sites <- rownames(rel_ref)

  # Pre-allocate results
  p_traits <- numeric(length(sites))
  p_abund  <- numeric(length(sites))
  ses_traits <- numeric(length(sites))
  ses_abund  <- numeric(length(sites))

  for (i in seq_along(sites)) {
    sitio <- sites[i]
    c_ref <- as.numeric(rel_ref[sitio, ])
    c_comp <- as.numeric(rel_comp[sitio, ])
    obs_dist <- sqrt(sum((x$CWM_ref[sitio, ] - x$CWM_comp[sitio, ])^2))

    # MODEL 1: Trait Permutation (Neutral Drift)
    sim_traits <- replicate(n_perm, {
      F_perm <- F_mat[sample(nrow(F_mat)), ]
      C1 <- c_ref %*% F_perm
      C2 <- c_comp %*% F_perm
      sqrt(sum((C1 - C2)^2))
    })

    # MODEL 2: Abundance Permutation (Structural Resilience)
    sim_abund <- replicate(n_perm, {
      c_comp_perm <- sample(c_comp)
      C2_perm <- c_comp_perm %*% F_mat
      sqrt(sum((x$CWM_ref[sitio, ] - C2_perm)^2))
    })

    # P-values
    p_traits[i] <- (sum(sim_traits >= obs_dist) + 1) / (n_perm + 1)
    p_abund[i]  <- (sum(sim_abund >= obs_dist) + 1) / (n_perm + 1)

    # Standardized Effect Sizes (SES) con parche de división por cero
    sd_traits <- sd(sim_traits)
    sd_abund <- sd(sim_abund)

    ses_traits[i] <- if (is.na(sd_traits) || sd_traits == 0) 0 else (obs_dist - mean(sim_traits)) / sd_traits
    ses_abund[i]  <- if (is.na(sd_abund) || sd_abund == 0) 0 else (obs_dist - mean(sim_abund)) / sd_abund
  }

  x$null_models <- data.frame(
    Contrast = sites,
    Obs_DeltaC = x$rDelta_C * (x$D_max / 100),
    p_Traits = p_traits,
    SES_Traits = ses_traits,
    p_Abund = p_abund,
    SES_Abund = ses_abund,
    stringsAsFactors = FALSE
  )

  cat("Done! Null models attached to object.\n")
  return(x)
}

#' @export
asc_null.ascfcd_pw <- function(x, n_perm = 999, seed = 42) {
  if (!is.null(seed)) set.seed(seed)

  n_pairs <- nrow(x$pairwise_results)
  cat(sprintf("\nRunning %d permutations for %d pairwise links...\n", n_perm, n_pairs))

  F_mat <- x$F_matrix
  rel_abund <- x$rel_abund
  df_res <- x$pairwise_results

  p_traits <- numeric(n_pairs)
  p_abund  <- numeric(n_pairs)
  ses_traits <- numeric(n_pairs)
  ses_abund  <- numeric(n_pairs)

  for (i in 1:n_pairs) {
    s1 <- df_res$Community_A[i]
    s2 <- df_res$Community_B[i]

    c1 <- as.numeric(rel_abund[s1, ])
    c2 <- as.numeric(rel_abund[s2, ])
    obs_dist <- sqrt(sum((x$CWM_matrix[s1, ] - x$CWM_matrix[s2, ])^2))

    # Model 1
    sim_traits <- replicate(n_perm, {
      F_perm <- F_mat[sample(nrow(F_mat)), ]
      C1_perm <- c1 %*% F_perm
      C2_perm <- c2 %*% F_perm
      sqrt(sum((C1_perm - C2_perm)^2))
    })

    # Model 2
    sim_abund <- replicate(n_perm, {
      c2_perm <- sample(c2)
      C2_perm <- c2_perm %*% F_mat
      sqrt(sum((x$CWM_matrix[s1, ] - C2_perm)^2))
    })

    p_traits[i] <- (sum(sim_traits >= obs_dist) + 1) / (n_perm + 1)
    p_abund[i]  <- (sum(sim_abund >= obs_dist) + 1) / (n_perm + 1)

    # Standardized Effect Sizes (SES) con parche de división por cero
    sd_traits <- sd(sim_traits)
    sd_abund <- sd(sim_abund)

    ses_traits[i] <- if (is.na(sd_traits) || sd_traits == 0) 0 else (obs_dist - mean(sim_traits)) / sd_traits
    ses_abund[i]  <- if (is.na(sd_abund) || sd_abund == 0) 0 else (obs_dist - mean(sim_abund)) / sd_abund
  }

  x$pairwise_results$p_Traits <- p_traits
  x$pairwise_results$SES_Traits <- ses_traits
  x$pairwise_results$p_Abund <- p_abund
  x$pairwise_results$SES_Abund <- ses_abund

  cat("Done! Null models attached to pairwise results.\n")
  return(x)
}
