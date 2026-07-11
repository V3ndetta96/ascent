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
#' @note
#' **Normalization:** Delta FDis and Delta FRic are absolute differences, not
#' normalized by the regional pool. Values are not directly comparable across
#' studies with different species pools or trait scales. Only rDelta_C (position)
#' is normalized by Dmax_regional.
#'
#' **FRic = NA:** When a community has fewer species than retained PCoA axes,
#' the convex hull is geometrically undefined and FRic is set to \code{NA},
#' which propagates to Delta_FRic.
#'
#' **Temporal replicates:** When multiple rows share the same site and temporal
#' state (e.g., three reference plots), their relative abundances are averaged
#' (column means) before computing the centroid. This is equivalent to treating
#' replicates as a single pooled community. Intra-state variability is not
#' propagated to downstream metrics.
#'
#' @return An S3 object of class \code{ascfcd}.
#'
#' @examples
#' # Simulate deforestation impact on a bird community
#' traits <- data.frame(
#'   Mass  = c(15, 30, 60, 150, 400),
#'   Beak  = c(10, 15, 28,  45,  85),
#'   Diet  = factor(c(0, 1, 1, 1, 0))
#' )
#' rownames(traits) <- paste0("Sp", 1:5)
#'
#' abund <- rbind(
#'   Reference = c(0, 10, 25, 20, 5),
#'   Impacted  = c(30, 15,  5,  0, 0)
#' )
#'
#' res <- asc_paired(
#'   traits, abund,
#'   sites = c("Site1", "Site1"),
#'   time = c("Reference", "Impacted"),
#'   ref_time = "Reference"
#' )
#' summary(res)
#'
#' @export
asc_paired <- function(traits, abund, sites, time, ref_time, dist_method = "gower",
                       dim_retention = c("variance", "broken_stick"), var_tol = 0.80, na.rm = TRUE) {

  dim_retention <- match.arg(dim_retention)

  if(nrow(traits) != ncol(abund)) stop("Number of species in traits and abundances does not match.")
  abund_rel <- as.matrix(sweep(abund, 1, rowSums(abund, na.rm = na.rm), "/"))

  fs <- .build_functional_space(traits, dist_method, dim_retention, var_tol)
  F_mat <- fs$F_mat

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

    # Guard: zero displacement
    ss_vdir <- sum(v_dir^2)
    if (ss_vdir < 1e-12) {
      asc_contrib <- rep(0, length(v_dir))
      crit_axis <- NA_integer_
    } else {
      asc_contrib <- (v_dir^2 / ss_vdir) * 100
      crit_axis <- which.max(asc_contrib)
    }

    topo_ref <- .calc_core_topology(as.numeric(p_ref_mat), F_mat)
    topo_comp <- .calc_core_topology(as.numeric(p_comp_mat), F_mat)

    res_list[[s]] <- list(
      abs_dist = dist_abs,
      asc = asc_contrib,
      critical_axis = crit_axis,
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
    F_space = F_mat, k_retained = fs$k_retained,
    var_retained = fs$var_retained, axis_var = fs$axis_var,
    quality = fs$quality,
    original_traits = traits, original_abund = abund_rel
  )
  class(output) <- "ascfcd"
  return(output)
}
