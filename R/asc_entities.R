#' Identify Functional Entities (FEs) and Redundancy
#'
#' @description
#' Groups species into Functional Entities (FEs) based on their unique trait combinations.
#' It calculates the functional redundancy (number of species per FE) and maps species
#' to their corresponding FE.
#'
#' @param traits A data frame or matrix of species (rows) by functional traits (columns).
#'
#' @return An S3 object of class \code{asc_fe} containing:
#' \describe{
#'   \item{fe_summary}{A data frame with FE IDs, species richness (redundancy), and trait profiles.}
#'   \item{fe_species}{A list where each element contains the names of the species belonging to a specific FE.}
#'   \item{functional_redundancy}{Numeric vector with the number of species per FE.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' fe_info <- asc_entities(traits)
#' }
asc_entities <- function(traits) {

  traits_df <- as.data.frame(traits)

  # Crear una "firma" única para cada especie uniendo sus valores de rasgos
  signatures <- apply(traits_df, 1, paste, collapse = "_")
  unique_sigs <- unique(signatures)
  n_fe <- length(unique_sigs)

  # Generar nombres limpios (FE_1, FE_2, ...)
  fe_names <- paste0("FE_", seq_len(n_fe))
  names(fe_names) <- unique_sigs

  # Mapear las firmas a las especies
  fe_species <- list()
  for(i in seq_along(unique_sigs)) {
    sig <- unique_sigs[i]
    spp <- rownames(traits_df)[signatures == sig]
    fe_species[[fe_names[sig]]] <- spp
  }

  # Crear la tabla resumen extrayendo el perfil de rasgos
  first_spp_idx <- match(unique_sigs, signatures)
  fe_summary <- traits_df[first_spp_idx, , drop = FALSE]
  rownames(fe_summary) <- fe_names[unique_sigs]

  # Añadir la columna de redundancia (Riqueza de especies por FE)
  redundancy <- sapply(fe_species, length)
  fe_summary$Redundancy <- redundancy

  # Reordenar para que la Redundancia sea la primera columna
  fe_summary <- cbind(Redundancy = fe_summary$Redundancy,
                      fe_summary[, -ncol(fe_summary), drop = FALSE])

  res <- list(
    fe_summary = fe_summary,
    fe_species = fe_species,
    functional_redundancy = redundancy,
    n_entities = n_fe,
    n_species = nrow(traits_df)
  )

  class(res) <- c("asc_fe", "list")
  return(res)
}
