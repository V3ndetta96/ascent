#' Pairwise Functional Spatial Divergence
#'
#' @description
#' Calculates spatial functional divergence across three layers: position (Delta C), dispersion (FDis), and volume (FRic).
#'
#' @param traits A data frame or matrix of functional traits.
#' @param abund A data frame or matrix of species abundances.
#' @param dist_method Distance metric for the trait matrix. Default is "gower".
#' @param dim_retention Character. Method for dimensionality reduction.
#' @param var_tol Numeric. Proportion of cumulative variance to retain. Default is 0.80.
#' @param na.rm Logical. Remove missing values? Default is TRUE.
#'
#' @return An S3 object of class \code{ascfcd_pw}.
#' @export
asc_pairwise <- function(traits, abund, dist_method = "gower",
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
  Dmax_empirical <- max(stats::dist(cwm_global))
  if(Dmax_empirical == 0) Dmax_empirical <- 1

  communities <- rownames(abund_rel)
  combinations <- utils::combn(communities, 2, simplify = FALSE)

  results <- data.frame()
  directional_vectors <- list(); p_ref_list <- list(); p_comp_list <- list()

  for(i in seq_along(combinations)) {
    com_A <- combinations[[i]][1]
    com_B <- combinations[[i]][2]
    pair_name <- paste(com_A, com_B, sep = "_vs_")

    p_ref_mat <- abund_rel[com_A, , drop = FALSE]
    p_comp_mat <- abund_rel[com_B, , drop = FALSE]
    p_ref_list[[pair_name]] <- p_ref_mat; p_comp_list[[pair_name]] <- p_comp_mat

    cwm_A <- cwm_global[com_A, , drop = FALSE]
    cwm_B <- cwm_global[com_B, , drop = FALSE]

    v_dir <- cwm_B - cwm_A
    directional_vectors[[pair_name]] <- v_dir
    dist_abs <- sqrt(sum(v_dir^2))
    r_delta <- (dist_abs / Dmax_empirical) * 100

    topo_A <- calc_topology(as.numeric(p_ref_mat), F_mat)
    topo_B <- calc_topology(as.numeric(p_comp_mat), F_mat)

    res_row <- data.frame(
      Community_A = com_A, Community_B = com_B,
      Delta_C_abs = dist_abs, rDelta_C_pct = r_delta,
      Delta_FDis = topo_B$FDis - topo_A$FDis,
      Delta_FRic = topo_B$FRic - topo_A$FRic
    )
    results <- rbind(results, res_row)
  }

  output <- list(
    pairwise_results = results, Dmax_empirical = Dmax_empirical,
    directional_vectors = directional_vectors,
    p_ref = p_ref_list, p_comp = p_comp_list,
    cwm_global = cwm_global, F_space = F_mat, k_retained = k_valid,
    var_retained = total_var_retained, axis_var = axis_var,
    original_traits = traits, original_abund = abund_rel
  )
  class(output) <- "ascfcd_pw"
  return(output)
}
