#' Print summary for ascfcd paired objects
#'
#' @param object An object of class ascfcd
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
summary.ascfcd <- function(object, ...) {
  cat("\n======================================================\n")
  cat("   ASC-FCD: Functional Centroid Displacement Summary\n")
  cat("======================================================\n")
  cat(sprintf("Functional Space: %d traits collapsed to %d effective axes.\n",
              ncol(object$traits), length(object$Var_j)))
  cat(sprintf("Maximum Functional Distance (D_max): %.4f\n", object$D_max))
  cat("------------------------------------------------------\n")

  site_names <- names(object$rDelta_C)
  has_null <- !is.null(object$null_models)

  for(s in site_names){
    base_str <- sprintf("Paired Contrast: %-10s | rDelta_C: %5.2f%% | Crit. Axis: %d (%.1f%% of change)",
                        s, object$rDelta_C[s], object$critical_axes[s], object$ASC_j_list[[s]][object$critical_axes[s]])

    if(has_null) {
      row_null <- object$null_models[object$null_models$Contrast == s, ]
      cat(sprintf("%s\n    -> p-val (Traits): %.3f | p-val (Abund): %.3f\n",
                  base_str, row_null$p_Traits, row_null$p_Abund))
    } else {
      cat(base_str, "\n")
    }
  }
  cat("======================================================\n")
}

#' Print summary for ascfcd pairwise objects
#'
#' @param object An object of class ascfcd_pw
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
summary.ascfcd_pw <- function(object, ...) {
  cat("\n======================================================\n")
  cat("   ASC-FCD: Pairwise Functional Displacement Summary\n")
  cat("======================================================\n")
  cat(sprintf("Functional Space: %d traits collapsed to %d effective axes.\n",
              ncol(object$traits), length(object$Var_j)))
  cat(sprintf("Communities compared: %d | Total pairwise links: %d\n",
              nrow(object$CWM_matrix), nrow(object$pairwise_results)))
  cat("------------------------------------------------------\n")
  cat("Top 5 largest functional shifts:\n")

  df_sorted <- object$pairwise_results[order(-object$pairwise_results$rDelta_C_pct), ]
  top_n <- min(5, nrow(df_sorted))
  has_null <- "p_Traits" %in% colnames(df_sorted)

  for(i in 1:top_n){
    row <- df_sorted[i, ]
    base_str <- sprintf("%-10s vs %-10s | rDelta_C: %5.2f%% | Crit. Axis: %d (%.1f%%)",
                        row$Community_A, row$Community_B,
                        row$rDelta_C_pct, row$Critical_Axis, row$ASC_Critical_pct)

    if(has_null) {
      cat(sprintf("%s\n    -> p-val (Traits): %.3f | p-val (Abund): %.3f\n",
                  base_str, row$p_Traits, row$p_Abund))
    } else {
      cat(base_str, "\n")
    }
  }
  cat("======================================================\n")
}

#' Print summary for asc_trans objects
#'
#' @param object An object of class asc_trans
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
summary.asc_trans <- function(object, ...) {
  cat("\n======================================================\n")
  cat("     ASC-FCD: Functional Traits & Species Transitions\n")
  cat("======================================================\n")

  print_spp <- function(spp_list, label) {
    cat(sprintf("\n%s:\n", label))
    if (length(spp_list) > 0) {
      for (sp in names(spp_list)) {
        cat(sprintf("  -> %-15s | Abund. local: %5.1f%% -> %5.1f%%\n",
                    sp, spp_list[[sp]]$Antes * 100, spp_list[[sp]]$Despues * 100))
      }
    } else {
      cat("  (Ninguna especie local cumple esta condicion)\n")
    }
  }

  contrasts <- names(object)

  for (cnt in contrasts) {
    cat(sprintf("\n>>> CONTRASTE / SITIO: %s (Eje Critico: %d)\n", toupper(cnt), object[[cnt]]$critical_axis))
    cat("------------------------------------------------------\n")
    cat("DIRECCION DE LOS RASGOS LIDERES (Basado en CWM):\n")

    traits_det <- object[[cnt]]$traits_evaluated
    for (tr in names(traits_det)) {
      t_data <- traits_det[[tr]]
      verbo_t <- ifelse(t_data$Dif > 0, "AUMENTO   [+]", "DISMINUYO [-]")

      if(t_data$is_binary) {
        cat(sprintf("  %s %-25s | %5.1f%% -> %5.1f%% (Dif: %5.1f%%)\n",
                    verbo_t, tr, t_data$Antes * 100, t_data$Despues * 100, t_data$Dif * 100))
      } else {
        cat(sprintf("  %s %-25s | %5.2f -> %5.2f (Dif: %5.2f)\n",
                    verbo_t, tr, t_data$Antes, t_data$Despues, t_data$Dif))
      }
    }

    print_spp(object[[cnt]]$win_strict,  "[++] GANADORAS ESTRICTAS (Impulsan TODOS los rasgos principales)")
    print_spp(object[[cnt]]$win_partial, "[+] GANADORAS PARCIALES (Impulsan AL MENOS 1 rasgo principal)")
    print_spp(object[[cnt]]$lose_strict, "[--] PERDEDORAS ESTRICTAS (Retraen TODOS los rasgos principales)")
    print_spp(object[[cnt]]$lose_partial, "[-] PERDEDORAS PARCIALES (Retraen AL MENOS 1 rasgo principal)")

    cat("======================================================\n")
  }
}

#' Print summary for asc_fe objects
#'
#' @param object An object of class asc_fe
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
summary.asc_fe <- function(object, ...) {
  cat("\n======================================================\n")
  cat("       ASC-FCD: Functional Entities Summary\n")
  cat("======================================================\n")
  cat(sprintf("Total Species Analyzed: %d\n", object$n_species))
  cat(sprintf("Unique Functional Entities (FEs): %d\n", object$n_entities))
  cat(sprintf("Average Redundancy: %.2f species per FE\n", mean(object$functional_redundancy)))
  cat("------------------------------------------------------\n")
  cat("Top 5 Most Redundant FEs, Species & Traits:\n")

  # Ordenar de mayor a menor redundancia
  df_sorted <- object$fe_summary[order(-object$fe_summary$Redundancy), ]
  top_n <- min(5, nrow(df_sorted))

  # Extraer la matriz de rasgos pura (sin la columna Redundancy)
  traits_df <- df_sorted[, -1, drop = FALSE]

  # PASO CLAVE: Identificar automáticamente cuáles rasgos son estrictamente binarios
  is_binary <- sapply(traits_df, function(col) {
    is.numeric(col) && all(col %in% c(0, 1, NA))
  })

  for(i in 1:top_n) {
    fe_id <- rownames(df_sorted)[i]
    redun <- df_sorted$Redundancy[i]

    # Extraer las especies
    spp_list <- paste(object$fe_species[[fe_id]], collapse = ", ")

    # Extraer y formatear los rasgos inteligentemente
    fe_vals <- traits_df[i, , drop = FALSE]
    rasgos_str_list <- c()

    for(j in seq_along(fe_vals)) {
      t_name <- colnames(fe_vals)[j]
      t_val <- fe_vals[1, j]

      if(is_binary[j]) {
        # LÓGICA CLÁSICA: Si es binario y vale 1, solo mostramos el nombre
        if(!is.na(t_val) && t_val == 1) {
          rasgos_str_list <- c(rasgos_str_list, t_name)
        }
      } else {
        # LÓGICA NUEVA: Si es continuo/categórico, mostramos "Nombre: Valor"
        if(is.numeric(t_val)) {
          rasgos_str_list <- c(rasgos_str_list, sprintf("%s: %.2f", t_name, t_val))
        } else {
          rasgos_str_list <- c(rasgos_str_list, sprintf("%s: %s", t_name, as.character(t_val)))
        }
      }
    }

    if (length(rasgos_str_list) > 0) {
      # Usamos una barra vertical para separar los rasgos continuos limpiamente
      rasgos_str <- paste(rasgos_str_list, collapse = " | ")
    } else {
      rasgos_str <- "(No se detectaron rasgos distintivos)"
    }

    cat(sprintf("\n  [%s] Redundancy: %d species\n", fe_id, redun))
    cat(sprintf("   -> Species: %s\n", spp_list))
    cat(sprintf("   -> Traits : %s\n", rasgos_str))
  }

  cat("\n======================================================\n")
  cat("* Tip: To see the species for ALL entities, type `your_object$fe_species`\n")
  cat("* Tip: To see the full trait matrix, type `your_object$fe_summary`\n")
  cat("======================================================\n")
}
