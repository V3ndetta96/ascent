#' Baseline Functional Topology of Entities
#'
#' @description
#' Calculates the absolute functional centroid (CWM), functional dispersion (FDis),
#' and convex hull volume (FRic) for each community or site independently.
#'
#' @param traits A data frame or matrix of functional traits.
#' @param abund A data frame or matrix of species abundances.
#' @param dist_method Distance metric for the trait matrix. Default is "gower".
#' @param dim_retention Character. Method for dimensionality reduction. Options are \code{"variance"} or \code{"broken_stick"}.
#' @param var_tol Numeric. Proportion of cumulative variance to retain. Default is 0.80.
#' @param na.rm Logical. Remove missing values? Default is TRUE.
#'
#' @return An S3 object of class \code{ascfcd_base}.
#' @export
asc_baseline <- function(traits, abund, dist_method = "gower",
                         dim_retention = c("variance", "broken_stick"), var_tol = 0.80, na.rm = TRUE) {

  dim_retention <- match.arg(dim_retention)
  if (!requireNamespace("geometry", quietly = TRUE)) {
    stop("The 'geometry' package is required for FRic calculations.")
  }

  if(nrow(traits) != ncol(abund)) stop("Number of species in traits and abundances does not match.")
  abund_rel <- as.matrix(sweep(abund, 1, rowSums(abund, na.rm = na.rm), "/"))

  d_mat <- cluster::daisy(traits, metric = dist_method)
  pcoa_res <- suppressWarnings(stats::cmdscale(d_mat, k = nrow(traits) - 1, eig = TRUE))

  eig_raw <- pcoa_res$eig
  eig_pos <- eig_raw[eig_raw > 1e-8]
  rel_var <- eig_pos / sum(eig_pos)
  cum_var <- cumsum(rel_var)

  if (dim_retention == "variance") {
    k_valid <- min(which(cum_var >= var_tol))
  } else if (dim_retention == "broken_stick") {
    n_eig <- length(eig_pos)
    bs_expected <- sapply(1:n_eig, function(k) sum(1 / (k:n_eig)) / n_eig)
    k_valid <- sum(rel_var > bs_expected)
  }
  if(k_valid < 2) k_valid <- 2

  axis_var <- rel_var[1:k_valid]
  total_var_retained <- cum_var[k_valid]
  F_mat <- pcoa_res$points[, 1:k_valid, drop = FALSE]
  colnames(F_mat) <- paste0("Axis_", 1:k_valid)

  # Motor Topológico Base
  calc_topo <- function(p_vec, f_space) {
    idx <- p_vec > 0
    if(sum(idx) < 2) return(list(FDis = 0, FRic = 0))
    p_pres <- p_vec[idx] / sum(p_vec[idx])
    f_pres <- f_space[idx, , drop = FALSE]
    cent <- colSums(p_pres * f_pres)
    fdis <- sum(p_pres * sqrt(rowSums(sweep(f_pres, 2, cent, "-")^2)))

    fric <- 0
    if (sum(idx) > ncol(f_space)) {
      fric <- tryCatch({ geometry::convhulln(f_pres, options = "FA")$vol }, error = function(e) 0)
    }
    return(list(FDis = fdis, FRic = fric))
  }

  cwm_global <- abund_rel %*% F_mat
  entities <- rownames(abund_rel)

  res_df <- data.frame(Entity = entities)
  for(j in 1:k_valid) {
    res_df[[paste0("CWM_Axis_", j)]] <- cwm_global[, j]
  }

  fdis_vals <- numeric(length(entities))
  fric_vals <- numeric(length(entities))

  for(i in seq_along(entities)) {
    topo <- calc_topo(abund_rel[i, ], F_mat)
    fdis_vals[i] <- topo$FDis
    fric_vals[i] <- topo$FRic
  }

  res_df$FDis <- fdis_vals
  res_df$FRic <- fric_vals

  output <- list(
    entities_results = res_df,
    cwm_global = cwm_global,
    F_space = F_mat,
    k_retained = k_valid,
    var_retained = total_var_retained,
    axis_var = axis_var,
    original_traits = traits,
    original_abund = abund_rel
  )
  class(output) <- "ascfcd_base"
  return(output)
}
