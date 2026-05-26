##' Compute Pairwise ASC-FCD Network
#'
#' @description
#' Computes the functional displacement for all possible pairwise combinations
#' of communities (rows) in the abundance matrix.
#'
#' @param traits A data frame or matrix of species traits (rows = species, columns = traits).
#' @param abund A data frame or matrix of species abundances (rows = communities, columns = species).
#' @param dist_method Distance method for traits (default: "jaccard").
#'
#' @return An object of class \code{ascfcd_pw}.
#' @export
#'
#' @importFrom vegan vegdist
#' @importFrom stats dist cmdscale
asc_pairwise <- function(traits, abund, dist_method = "jaccard") {

  # 1. Alinear matrices
  spp_abund <- colnames(abund)
  spp_traits <- rownames(traits)
  common_spp <- intersect(spp_abund, spp_traits)

  if(length(common_spp) == 0) stop("No common species between traits and abundance matrices.")

  abund <- abund[, common_spp, drop = FALSE]
  traits <- traits[common_spp, , drop = FALSE]

  # 2. Relativizar abundancias
  rel_abund <- sweep(abund, 1, rowSums(abund), "/")
  rel_abund[is.na(rel_abund)] <- 0

  # 3. PCoA (Con la corrección robusta de dimensiones)
  dist_func <- vegan::vegdist(traits, method = dist_method)
  pcoa_res <- stats::cmdscale(dist_func, k = ncol(traits) - 1, eig = TRUE)

  k_returned <- ncol(pcoa_res$points)
  eig_returned <- pcoa_res$eig[1:k_returned]

  valid_axes <- eig_returned > 1e-8
  if(sum(valid_axes) == 0) stop("No positive eigenvalues in PCoA.")

  F_matrix <- pcoa_res$points[, valid_axes, drop = FALSE]
  Var_j <- eig_returned[valid_axes] / sum(pcoa_res$eig[pcoa_res$eig > 0])

  # 4. Geometría espacial
  D_max <- sqrt(sum((apply(F_matrix, 2, max) - apply(F_matrix, 2, min))^2))
  CWM_matrix <- as.matrix(rel_abund) %*% F_matrix

  n_sites <- nrow(CWM_matrix)
  site_names <- rownames(CWM_matrix)

  # 5. Calcular todos los enlaces posibles (Combinatoria pura)
  res_list <- list()
  counter <- 1

  for(i in 1:(n_sites - 1)) {
    for(j in (i + 1):n_sites) {
      s1 <- site_names[i]
      s2 <- site_names[j]

      delta_j <- CWM_matrix[s2, ] - CWM_matrix[s1, ]
      Delta_C <- sqrt(sum(delta_j^2))
      rDelta_C <- (Delta_C / D_max) * 100

      d_j <- abs(delta_j)
      asc_j <- (d_j * Var_j) / sum(d_j * Var_j) * 100
      crit_axis <- which.max(asc_j)

      res_list[[counter]] <- data.frame(
        Community_A = s1,
        Community_B = s2,
        Delta_C_abs = Delta_C,
        rDelta_C_pct = rDelta_C,
        Critical_Axis = crit_axis,
        ASC_Critical_pct = asc_j[crit_axis],
        stringsAsFactors = FALSE
      )
      counter <- counter + 1
    }
  }

  df_results <- do.call(rbind, res_list)

  res <- list(
    traits = traits,
    rel_abund = rel_abund,
    CWM_matrix = CWM_matrix,
    F_matrix = F_matrix,
    Var_j = Var_j,
    D_max = D_max,
    pairwise_results = df_results
  )

  class(res) <- c("ascfcd_pw", "list")
  return(res)
}
