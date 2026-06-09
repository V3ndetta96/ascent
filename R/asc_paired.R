#' Paired Functional Multidimensional Restructuring
#'
#' @description
#' Calculates centroid displacement (Layer 1), dispersion shift (Layer 2 - FDis),
#' and volume shift (Layer 3 - FRic) between paired states.
#'
#' @param traits A data frame or matrix of functional traits.
#' @param abund A data frame or matrix of species abundances.
#' @param sites A vector indicating the site identifier.
#' @param time A vector indicating the temporal state.
#' @param ref_time A character string specifying the reference state.
#' @param dist_method Distance metric for the trait matrix. Default is "gower".
#' @param dim_retention Character. Method for dimensionality reduction.
#' @param var_tol Numeric. Proportion of cumulative variance to retain. Default is 0.80.
#' @param na.rm Logical. Remove missing values? Default is TRUE.
#'
#' @return An S3 object of class \code{ascfcd}.
#' @export
asc_paired <- function(traits, abund, sites, time, ref_time, dist_method = "gower",
                       dim_retention = c("variance", "broken_stick"), var_tol = 0.80, na.rm = TRUE) {

  dim_retention <- match.arg(dim_retention)
  if (!requireNamespace("geometry", quietly = TRUE)) {
    stop("The 'geometry' package is required for FRic calculations. Please install it.")
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

  # Motor Topológico Multicapa (FDis y FRic)
  calc_topology <- function(abund_vec, F_space) {
    idx <- abund_vec > 0
    if(sum(idx) < 2) return(list(FDis = 0, FRic = 0))

    p_pres <- abund_vec[idx] / sum(abund_vec[idx])
    F_pres <- F_space[idx, , drop = FALSE]

    centroid <- colSums(p_pres * F_pres)
    dist_to_c <- sqrt(rowSums(sweep(F_pres, 2, centroid, "-")^2))
    fdis <- sum(p_pres * dist_to_c)

    fric <- 0
    if (sum(idx) > ncol(F_space)) {
      fric <- tryCatch({
        geometry::convhulln(F_pres, options = "FA")$vol
      }, error = function(e) 0)
    }
    return(list(FDis = fdis, FRic = fric))
  }

  cwm_global <- abund_rel %*% F_mat
  Dmax_regional <- max(stats::dist(F_mat))
  if(Dmax_regional == 0) Dmax_regional <- 1

  unique_sites <- unique(sites)
  res_list <- list()
  rDelta_C <- numeric(length(unique_sites))
  names(rDelta_C) <- unique_sites
  directional_vectors <- list(); p_ref_list <- list(); p_comp_list <- list()
  cwm_ref_list <- list(); cwm_comp_list <- list()

  for (s in unique_sites) {
    idx <- which(sites == s)
    abund_site <- abund_rel[idx, , drop = FALSE]
    time_site <- time[idx]

    idx_ref <- which(time_site == ref_time)
    idx_comp <- which(time_site != ref_time)

    if(length(idx_ref) == 0 || length(idx_comp) == 0) { rDelta_C[s] <- NA; next }

    p_ref_mat <- if(nrow(abund_site[idx_ref, , drop=FALSE]) > 1) t(as.matrix(colMeans(abund_site[idx_ref, , drop=FALSE]))) else abund_site[idx_ref, , drop=FALSE]
    p_comp_mat <- if(nrow(abund_site[idx_comp, , drop=FALSE]) > 1) t(as.matrix(colMeans(abund_site[idx_comp, , drop=FALSE]))) else abund_site[idx_comp, , drop=FALSE]

    p_ref_list[[s]] <- p_ref_mat; p_comp_list[[s]] <- p_comp_mat
    cwm_ref <- p_ref_mat %*% F_mat; cwm_comp <- p_comp_mat %*% F_mat
    cwm_ref_list[[s]] <- cwm_ref; cwm_comp_list[[s]] <- cwm_comp

    v_dir <- cwm_comp - cwm_ref
    directional_vectors[[s]] <- v_dir
    dist_abs <- sqrt(sum(v_dir^2))
    rDelta_C[s] <- (dist_abs / Dmax_regional) * 100

    topo_ref <- calc_topology(as.numeric(p_ref_mat), F_mat)
    topo_comp <- calc_topology(as.numeric(p_comp_mat), F_mat)

    res_list[[s]] <- list(
      abs_dist = dist_abs,
      asc = (v_dir^2 / sum(v_dir^2)) * 100,
      critical_axis = which.max((v_dir^2 / sum(v_dir^2)) * 100),
      Delta_FDis = topo_comp$FDis - topo_ref$FDis,
      Delta_FRic = topo_comp$FRic - topo_ref$FRic
    )
  }

  valid_sites <- !is.na(rDelta_C)
  output <- list(
    rDelta_C = rDelta_C[valid_sites],
    Dmax_regional = Dmax_regional,
    directional_vectors = directional_vectors,
    p_ref = p_ref_list, p_comp = p_comp_list,
    cwm_ref = cwm_ref_list, cwm_comp = cwm_comp_list,
    site_results = res_list[valid_sites],
    F_space = F_mat, k_retained = k_valid,
    var_retained = total_var_retained, axis_var = axis_var,
    original_traits = traits, original_abund = abund_rel
  )
  class(output) <- "ascfcd"
  return(output)
}
