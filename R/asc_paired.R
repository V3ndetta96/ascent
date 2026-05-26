#' Compute Paired ASC-FCD
#'
#' @description
#' Evaluates functional centroid shifts between paired communities (e.g., Before/After,
#' Control/Impact) using a single abundance matrix and grouping vectors.
#'
#' @param traits A data frame or matrix of species traits (rows = species, columns = traits).
#' @param abund A data frame or matrix of species abundances (rows = samples, columns = species).
#' @param sites A character vector indicating the site/plot identity for each row in \code{abund}.
#' @param time A character vector indicating the temporal or experimental group for each row. Must have exactly two levels.
#' @param ref_time Character. The specific value in \code{time} that acts as the baseline/reference. If NULL, uses the first unique value.
#' @param dist_method Distance method for traits (default: "jaccard").
#'
#' @return An object of class \code{ascfcd}.
#' @export
#'
#' @importFrom vegan vegdist
#' @importFrom stats dist cmdscale
asc_paired <- function(traits, abund, sites, time, ref_time = NULL, dist_method = "jaccard") {

  # 1. Validaciones
  if (nrow(abund) != length(sites) || nrow(abund) != length(time)) {
    stop("The length of 'sites' and 'time' must match the number of rows in 'abund'.")
  }

  time_levels <- unique(time)
  if (length(time_levels) != 2) {
    stop(sprintf("The 'time' vector must have exactly two unique conditions. Found: %d", length(time_levels)))
  }

  if (is.null(ref_time)) {
    ref_time <- time_levels[1]
    comp_time <- time_levels[2]
    message(sprintf("No ref_time provided. Using '%s' as Reference and '%s' as Comparison.", ref_time, comp_time))
  } else {
    if (!(ref_time %in% time_levels)) {
      stop(sprintf("The ref_time '%s' is not present in the 'time' vector.", ref_time))
    }
    comp_time <- setdiff(time_levels, ref_time)
  }

  # 2. Separar y Alinear Matrices automáticamente
  abund_ref <- abund[time == ref_time, , drop = FALSE]
  sites_ref <- sites[time == ref_time]

  abund_comp <- abund[time == comp_time, , drop = FALSE]
  sites_comp <- sites[time == comp_time]

  common_sites <- intersect(sites_ref, sites_comp)
  if (length(common_sites) == 0) {
    stop("No matching sites found between the two time periods.")
  }

  rownames(abund_ref) <- sites_ref
  rownames(abund_comp) <- sites_comp

  abund_ref <- abund_ref[common_sites, , drop = FALSE]
  abund_comp <- abund_comp[common_sites, , drop = FALSE]

  # 3. Alinear con los rasgos
  spp_abund <- colnames(abund_ref)
  spp_traits <- rownames(traits)
  common_spp <- intersect(spp_abund, spp_traits)

  if(length(common_spp) == 0) stop("No common species between traits and abundance matrices.")

  abund_ref <- abund_ref[, common_spp, drop = FALSE]
  abund_comp <- abund_comp[, common_spp, drop = FALSE]
  traits <- traits[common_spp, , drop = FALSE]

  # 4. Abundancias Relativas
  rel_ref <- sweep(abund_ref, 1, rowSums(abund_ref), "/")
  rel_ref[is.na(rel_ref)] <- 0

  rel_comp <- sweep(abund_comp, 1, rowSums(abund_comp), "/")
  rel_comp[is.na(rel_comp)] <- 0

  # 5. Cálculo del PCoA
  dist_func <- vegan::vegdist(traits, method = dist_method)
  pcoa_res <- stats::cmdscale(dist_func, k = ncol(traits) - 1, eig = TRUE)

  # Ajustar el filtro a las dimensiones reales generadas
  k_returned <- ncol(pcoa_res$points)
  eig_returned <- pcoa_res$eig[1:k_returned]

  valid_axes <- eig_returned > 1e-8
  if(sum(valid_axes) == 0) stop("No positive eigenvalues in PCoA.")

  F_matrix <- pcoa_res$points[, valid_axes, drop = FALSE]

  # Calcular la varianza respecto a la inercia positiva total del sistema
  Var_j <- eig_returned[valid_axes] / sum(pcoa_res$eig[pcoa_res$eig > 0])

  # 6. Distancia Máxima y Desplazamientos
  D_max <- sqrt(sum((apply(F_matrix, 2, max) - apply(F_matrix, 2, min))^2))

  CWM_ref <- as.matrix(rel_ref) %*% F_matrix
  CWM_comp <- as.matrix(rel_comp) %*% F_matrix

  delta_CWM <- CWM_comp - CWM_ref
  Delta_C <- sqrt(rowSums(delta_CWM^2))
  rDelta_C <- (Delta_C / D_max) * 100
  names(rDelta_C) <- common_sites

  ASC_j_list <- list()
  critical_axes <- numeric(length(common_sites))
  names(critical_axes) <- common_sites

  for(i in seq_along(common_sites)) {
    sitio <- common_sites[i]
    delta_j <- delta_CWM[sitio, ]
    d_j <- abs(delta_j)

    asc_j <- (d_j * Var_j) / sum(d_j * Var_j) * 100
    ASC_j_list[[sitio]] <- asc_j
    critical_axes[sitio] <- which.max(asc_j)
  }

  res <- list(
    traits = traits,
    rel_ref = rel_ref,
    rel_comp = rel_comp,
    CWM_ref = CWM_ref,
    CWM_comp = CWM_comp,
    F_matrix = F_matrix,
    Var_j = Var_j,
    D_max = D_max,
    rDelta_C = rDelta_C,
    ASC_j_list = ASC_j_list,
    critical_axes = critical_axes,
    null_models = NULL
  )

  class(res) <- c("ascfcd", "list")
  return(res)
}
