#' Identify Functional Entity Transitions and Driver Species
#'
#' @description
#' Decomposes the critical axes of change into specific trait shifts and identifies
#' the taxonomic drivers using Community Weighted Means (CWM) and Pearson correlation.
#' Automatically handles binary and continuous traits.
#'
#' @param x An object of class `ascfcd` (paired) or `ascfcd_pw` (pairwise).
#' @param top_n Integer. Number of top contributing traits to consider (default: 3).
#'
#' @return A structured list of class `asc_trans`.
#' @export
asc_transitions <- function(x, top_n = 3) {
  UseMethod("asc_transitions")
}

#' @export
asc_transitions.ascfcd <- function(x, top_n = 3) {
  cat("\nExtracting functional transitions for paired communities...\n")

  F_matrix <- x$F_matrix
  matriz_rasgos <- as.data.frame(x$traits)
  sites <- names(x$rDelta_C)

  # Quedarnos solo con las columnas numéricas para el CWM y Correlación
  is_num <- sapply(matriz_rasgos, is.numeric)
  mat_num <- matriz_rasgos[, is_num, drop = FALSE]

  results_list <- list()

  for (sitio in sites) {
    eje_critico <- x$critical_axes[sitio]

    abund_antes <- as.numeric(x$rel_ref[sitio, ])
    names(abund_antes) <- colnames(x$rel_ref)
    abund_desp  <- as.numeric(x$rel_comp[sitio, ])
    names(abund_desp) <- colnames(x$rel_comp)

    coords_eje <- F_matrix[, eje_critico]

    # 1. Identificar la importancia del rasgo usando Correlación de Pearson
    correlaciones <- numeric(ncol(mat_num))
    names(correlaciones) <- colnames(mat_num)
    for(i in 1:ncol(mat_num)) {
      if(sd(mat_num[, i], na.rm = TRUE) > 0) {
        correlaciones[i] <- abs(stats::cor(mat_num[, i], coords_eje, use = "pairwise.complete.obs"))
      } else {
        correlaciones[i] <- 0
      }
    }

    top_rasgos <- names(sort(correlaciones, decreasing = TRUE)[1:top_n])

    rasgos_aum <- c()
    rasgos_dis <- c()
    traits_summary <- list()

    # Matriz indicadora temporal para agrupar ganadoras/perdedoras
    mat_ind <- matrix(0, nrow = nrow(mat_num), ncol = length(top_rasgos))
    rownames(mat_ind) <- rownames(mat_num)
    colnames(mat_ind) <- top_rasgos

    for (rasgo in top_rasgos) {
      rasgo_vals <- mat_num[, rasgo]
      is_binary <- all(rasgo_vals %in% c(0, 1, NA))

      # 2. Calcular la Media Ponderada de la Comunidad (CWM)
      cwm_antes <- sum(abund_antes * rasgo_vals, na.rm = TRUE)
      cwm_desp  <- sum(abund_desp * rasgo_vals, na.rm = TRUE)
      cambio <- cwm_desp - cwm_antes

      # 3. Identificar quién "posee" el rasgo conductor
      if(is_binary) {
        spp_con_rasgo <- rownames(mat_num)[rasgo_vals == 1 & !is.na(rasgo_vals)]
      } else {
        spp_con_rasgo <- rownames(mat_num)[rasgo_vals > mean(rasgo_vals, na.rm = TRUE) & !is.na(rasgo_vals)]
      }

      mat_ind[spp_con_rasgo, rasgo] <- 1

      traits_summary[[rasgo]] <- data.frame(Antes = cwm_antes, Despues = cwm_desp,
                                            Dif = cambio, is_binary = is_binary)

      if (cambio > 0) { rasgos_aum <- c(rasgos_aum, rasgo) }
      else { rasgos_dis <- c(rasgos_dis, rasgo) }
    }

    spp_win_strict <- list(); spp_win_partial <- list()
    spp_lose_strict <- list(); spp_lose_partial <- list()

    if (length(rasgos_aum) > 0) {
      mat_aum <- mat_ind[, rasgos_aum, drop = FALSE]
      sums_aum <- rowSums(mat_aum)
      for (sp in rownames(mat_aum)) {
        if (abund_desp[sp] > abund_antes[sp]) {
          if (sums_aum[sp] == length(rasgos_aum)) {
            spp_win_strict[[sp]] <- data.frame(Antes = abund_antes[sp], Despues = abund_desp[sp])
          } else if (sums_aum[sp] >= 1) {
            spp_win_partial[[sp]] <- data.frame(Antes = abund_antes[sp], Despues = abund_desp[sp])
          }
        }
      }
    }

    if (length(rasgos_dis) > 0) {
      mat_dis <- mat_ind[, rasgos_dis, drop = FALSE]
      sums_dis <- rowSums(mat_dis)
      for (sp in rownames(mat_dis)) {
        if (abund_desp[sp] < abund_antes[sp]) {
          if (sums_dis[sp] == length(rasgos_dis)) {
            spp_lose_strict[[sp]] <- data.frame(Antes = abund_antes[sp], Despues = abund_desp[sp])
          } else if (sums_dis[sp] >= 1) {
            spp_lose_partial[[sp]] <- data.frame(Antes = abund_antes[sp], Despues = abund_desp[sp])
          }
        }
      }
    }

    results_list[[sitio]] <- list(
      critical_axis = eje_critico,
      traits_evaluated = traits_summary,
      win_strict = spp_win_strict,
      win_partial = spp_win_partial,
      lose_strict = spp_lose_strict,
      lose_partial = spp_lose_partial
    )
  }

  class(results_list) <- c("asc_trans", "list")
  return(results_list)
}

#' @export
asc_transitions.ascfcd_pw <- function(x, top_n = 3) {
  cat("\nExtracting functional transitions for pairwise network links...\n")

  F_matrix <- x$F_matrix
  matriz_rasgos <- as.data.frame(x$traits)
  df_res <- x$pairwise_results

  is_num <- sapply(matriz_rasgos, is.numeric)
  mat_num <- matriz_rasgos[, is_num, drop = FALSE]

  results_list <- list()

  for (i in 1:nrow(df_res)) {
    s1 <- df_res$Community_A[i]
    s2 <- df_res$Community_B[i]
    link_name <- paste(s1, "vs", s2)

    eje_critico <- df_res$Critical_Axis[i]

    abund_antes <- as.numeric(x$rel_abund[s1, ])
    names(abund_antes) <- colnames(x$rel_abund)
    abund_desp  <- as.numeric(x$rel_abund[s2, ])
    names(abund_desp) <- colnames(x$rel_abund)

    coords_eje <- F_matrix[, eje_critico]

    correlaciones <- numeric(ncol(mat_num))
    names(correlaciones) <- colnames(mat_num)
    for(j in 1:ncol(mat_num)) {
      if(sd(mat_num[, j], na.rm = TRUE) > 0) {
        correlaciones[j] <- abs(stats::cor(mat_num[, j], coords_eje, use = "pairwise.complete.obs"))
      } else {
        correlaciones[j] <- 0
      }
    }

    top_rasgos <- names(sort(correlaciones, decreasing = TRUE)[1:top_n])

    rasgos_aum <- c()
    rasgos_dis <- c()
    traits_summary <- list()

    mat_ind <- matrix(0, nrow = nrow(mat_num), ncol = length(top_rasgos))
    rownames(mat_ind) <- rownames(mat_num)
    colnames(mat_ind) <- top_rasgos

    for (rasgo in top_rasgos) {
      rasgo_vals <- mat_num[, rasgo]
      is_binary <- all(rasgo_vals %in% c(0, 1, NA))

      cwm_antes <- sum(abund_antes * rasgo_vals, na.rm = TRUE)
      cwm_desp  <- sum(abund_desp * rasgo_vals, na.rm = TRUE)
      cambio <- cwm_desp - cwm_antes

      if(is_binary) {
        spp_con_rasgo <- rownames(mat_num)[rasgo_vals == 1 & !is.na(rasgo_vals)]
      } else {
        spp_con_rasgo <- rownames(mat_num)[rasgo_vals > mean(rasgo_vals, na.rm = TRUE) & !is.na(rasgo_vals)]
      }

      mat_ind[spp_con_rasgo, rasgo] <- 1

      traits_summary[[rasgo]] <- data.frame(Antes = cwm_antes, Despues = cwm_desp,
                                            Dif = cambio, is_binary = is_binary)

      if (cambio > 0) { rasgos_aum <- c(rasgos_aum, rasgo) }
      else { rasgos_dis <- c(rasgos_dis, rasgo) }
    }

    spp_win_strict <- list(); spp_win_partial <- list()
    spp_lose_strict <- list(); spp_lose_partial <- list()

    if (length(rasgos_aum) > 0) {
      mat_aum <- mat_ind[, rasgos_aum, drop = FALSE]
      sums_aum <- rowSums(mat_aum)
      for (sp in rownames(mat_aum)) {
        if (abund_desp[sp] > abund_antes[sp]) {
          if (sums_aum[sp] == length(rasgos_aum)) {
            spp_win_strict[[sp]] <- data.frame(Antes = abund_antes[sp], Despues = abund_desp[sp])
          } else if (sums_aum[sp] >= 1) {
            spp_win_partial[[sp]] <- data.frame(Antes = abund_antes[sp], Despues = abund_desp[sp])
          }
        }
      }
    }

    if (length(rasgos_dis) > 0) {
      mat_dis <- mat_ind[, rasgos_dis, drop = FALSE]
      sums_dis <- rowSums(mat_dis)
      for (sp in rownames(mat_dis)) {
        if (abund_desp[sp] < abund_antes[sp]) {
          if (sums_dis[sp] == length(rasgos_dis)) {
            spp_lose_strict[[sp]] <- data.frame(Antes = abund_antes[sp], Despues = abund_desp[sp])
          } else if (sums_dis[sp] >= 1) {
            spp_lose_partial[[sp]] <- data.frame(Antes = abund_antes[sp], Despues = abund_desp[sp])
          }
        }
      }
    }

    results_list[[link_name]] <- list(
      critical_axis = eje_critico,
      traits_evaluated = traits_summary,
      win_strict = spp_win_strict,
      win_partial = spp_win_partial,
      lose_strict = spp_lose_strict,
      lose_partial = spp_lose_partial
    )
  }

  class(results_list) <- c("asc_trans", "list")
  return(results_list)
}
