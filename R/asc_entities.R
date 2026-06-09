#' Identify Functional Entities (Species Clustering)
#'
#' @description
#' Clusters species into discrete functional entities based on shared
#' morphological, physiological, or ecological traits.
#'
#' @param traits A data frame or matrix of functional traits.
#' @param dist_method Character. Distance metric. Default is \code{"gower"} to handle mixed data types.
#' @param hclust_method Character. Agglomeration method for hierarchical clustering. Default is \code{"ward.D2"}.
#' @param k Integer. Desired number of functional entities (clusters). If \code{NULL}, \code{h} must be provided.
#' @param h Numeric. Height at which to cut the dendrogram. If both \code{k} and \code{h} are \code{NULL}, defaults to \code{k = 3}.
#'
#' @return An S3 object of class \code{ascfcd_entities}.
#' @export
asc_entities <- function(traits, dist_method = "gower", hclust_method = "ward.D2", k = NULL, h = NULL) {

  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("The 'cluster' package is required for distance calculations.")
  }

  # 1. Cálculo de Distancias
  d_mat <- cluster::daisy(traits, metric = dist_method)

  # 2. Clustering Jerárquico
  hc <- stats::hclust(d_mat, method = hclust_method)

  # 3. Determinación del corte (agrupamiento)
  if (is.null(k) && is.null(h)) {
    warning("Neither 'k' (number of clusters) nor 'h' (cut height) were provided. Defaulting to k = 3 functional entities.")
    k <- 3
  }

  clusters <- stats::cutree(hc, k = k, h = h)

  # 4. Ensamblaje de resultados
  res_df <- data.frame(
    Species = rownames(traits),
    Entity_ID = paste0("Entity_", clusters),
    stringsAsFactors = FALSE
  )

  # Añadimos los rasgos originales para facilitar la inspección del usuario
  res_df <- cbind(res_df, traits)

  output <- list(
    entity_classification = res_df,
    distance_matrix = d_mat,
    hclust_obj = hc,
    method = hclust_method,
    metric = dist_method
  )

  class(output) <- "ascfcd_entities"
  return(output)
}
