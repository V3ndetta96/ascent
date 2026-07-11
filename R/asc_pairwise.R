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
#' @note
#' **Leverage direction:** In pairwise mode, leverage is computed from
#' Community_A toward Community_B. Reversing the pair order reverses the
#' leverage signs. The direction is determined by alphabetical ordering of
#' community names.
#'
#' **Normalization:** Delta FDis and Delta FRic are absolute differences (see
#' \code{\link{asc_paired}} for details).
#'
#' @return An S3 object of class \code{ascfcd_pw}.
#'
#' @examples
#' traits <- data.frame(
#'   Mass = c(15, 30, 60, 150, 400),
#'   Beak = c(10, 15, 28,  45,  85)
#' )
#' rownames(traits) <- paste0("Sp", 1:5)
#'
#' abund <- rbind(
#'   Forest  = c(0,  5, 25, 20, 10),
#'   Field   = c(30, 20, 10,  0,  0),
#'   Wetland = c(5, 10, 15, 15,  5)
#' )
#'
#' res_pw <- asc_pairwise(traits, abund, dist_method = "euclidean")
#' summary(res_pw)
#'
#' @export
asc_pairwise <- function(traits, abund, dist_method = "gower",
                         dim_retention = c("variance", "broken_stick"), var_tol = 0.80, na.rm = TRUE) {

  dim_retention <- match.arg(dim_retention)

  if(nrow(traits) != ncol(abund)) stop("Number of species in traits and abundances does not match.")
  abund_rel <- as.matrix(sweep(abund, 1, rowSums(abund, na.rm = na.rm), "/"))

  fs <- .build_functional_space(traits, dist_method, dim_retention, var_tol)
  F_mat <- fs$F_mat

  cwm_global <- abund_rel %*% F_mat
  Dmax_regional <- max(stats::dist(F_mat))
  if(Dmax_regional == 0) Dmax_regional <- 1

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
    r_delta <- (dist_abs / Dmax_regional) * 100

    topo_A <- .calc_core_topology(as.numeric(p_ref_mat), F_mat)
    topo_B <- .calc_core_topology(as.numeric(p_comp_mat), F_mat)

    res_row <- data.frame(
      Community_A = com_A, Community_B = com_B,
      Delta_C_abs = dist_abs, rDelta_C_pct = r_delta,
      Delta_FDis = topo_B$FDis - topo_A$FDis,
      Delta_FRic = topo_B$FRic - topo_A$FRic
    )
    results <- rbind(results, res_row)
  }

  output <- list(
    pairwise_results = results, Dmax_regional = Dmax_regional,
    directional_vectors = directional_vectors,
    p_ref = p_ref_list, p_comp = p_comp_list,
    cwm_global = cwm_global, F_space = F_mat, k_retained = fs$k_retained,
    var_retained = fs$var_retained, axis_var = fs$axis_var,
    quality = fs$quality,
    original_traits = traits, original_abund = abund_rel
  )
  class(output) <- "ascfcd_pw"
  return(output)
}
