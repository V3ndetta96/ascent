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
#'
#' @examples
#' traits <- data.frame(
#'   Mass = c(15, 30, 60, 150, 400),
#'   Beak = c(10, 15, 28,  45,  85)
#' )
#' rownames(traits) <- paste0("Sp", 1:5)
#'
#' abund <- rbind(
#'   Forest = c(0,  5, 25, 20, 10),
#'   Field  = c(30, 20, 10,  0,  0)
#' )
#'
#' base <- asc_baseline(traits, abund, dist_method = "euclidean")
#' summary(base)
#'
#' @export
asc_baseline <- function(traits, abund, dist_method = "gower",
                         dim_retention = c("variance", "broken_stick"), var_tol = 0.80, na.rm = TRUE) {

  dim_retention <- match.arg(dim_retention)

  if(nrow(traits) != ncol(abund)) stop("Number of species in traits and abundances does not match.")
  abund_rel <- as.matrix(sweep(abund, 1, rowSums(abund, na.rm = na.rm), "/"))

  fs <- .build_functional_space(traits, dist_method, dim_retention, var_tol)
  F_mat <- fs$F_mat
  k_valid <- fs$k_retained

  cwm_global <- abund_rel %*% F_mat
  entities <- rownames(abund_rel)

  res_df <- data.frame(Entity = entities)
  for(j in 1:k_valid) {
    res_df[[paste0("CWM_Axis_", j)]] <- cwm_global[, j]
  }

  fdis_vals <- numeric(length(entities))
  fric_vals <- numeric(length(entities))

  for(i in seq_along(entities)) {
    topo <- .calc_core_topology(abund_rel[i, ], F_mat)
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
    var_retained = fs$var_retained,
    axis_var = fs$axis_var,
    quality = fs$quality,
    original_traits = traits,
    original_abund = abund_rel
  )
  class(output) <- "ascfcd_base"
  return(output)
}
