#' Extract Functional Transitions and Species Leverage
#'
#' @description
#' Calculates the multidimensional Functional Leverage of each species,
#' quantifying their direct contribution to the positional shift (Layer 1 - Delta C)
#' of the ecosystem.
#'
#' @param x An object of class \code{ascfcd} or \code{ascfcd_pw}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A list of data frames containing species leverage for each contrast.
#' @export
asc_transitions <- function(x, ...) {
  if (!inherits(x, c("ascfcd", "ascfcd_pw"))) {
    stop("Input must be of class 'ascfcd' or 'ascfcd_pw'.")
  }

  F_mat <- x$F_space
  entities <- names(x$directional_vectors)
  res_list <- list()

  for (s in entities) {
    v_dir <- x$directional_vectors[[s]]
    mag_v <- sqrt(sum(v_dir^2))

    # Extraemos la demografía
    p_r <- as.numeric(x$p_ref[[s]])
    p_c <- as.numeric(x$p_comp[[s]])
    delta_p <- p_c - p_r

    if (mag_v == 0) {
      # Si no hay desplazamiento, el leverage es nulo
      df_lev <- data.frame(Species = rownames(F_mat), Delta_p = delta_p, Projection = 0, Leverage = 0)
    } else {
      # Vector unitario direccional
      u_dir <- v_dir / mag_v

      # Ubicamos el centroide de origen según el tipo de objeto
      if (inherits(x, "ascfcd")) {
        c_base <- as.numeric(x$cwm_ref[[s]])
      } else {
        com_A <- strsplit(s, "_vs_")[[1]][1]
        c_base <- as.numeric(x$cwm_global[com_A, ])
      }

      # Proyección ortogonal de cada especie sobre el vector direccional
      projections <- numeric(nrow(F_mat))
      for (i in seq_len(nrow(F_mat))) {
        projections[i] <- sum((F_mat[i, ] - c_base) * u_dir)
      }

      # Cálculo del Leverage Multidimensional
      leverage <- delta_p * projections

      df_lev <- data.frame(
        Species = rownames(F_mat),
        Delta_p = round(delta_p, 4),
        Projection = round(projections, 4),
        Leverage = round(leverage, 4)
      )
    }

    # Ordenamos por impacto absoluto (los mayores drivers arriba)
    df_lev <- df_lev[order(abs(df_lev$Leverage), decreasing = TRUE), ]
    rownames(df_lev) <- NULL

    res_list[[s]] <- list(species_leverage = df_lev)
  }

  return(res_list)
}
