#' Summary for Paired Functional Centroid Displacement
#'
#' @description
#' Provides a concise, multi-layer overview of the paired ASC-FCD analysis,
#' including null models and top functional drivers.
#'
#' @param object An object of class \code{ascfcd}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns a data frame with the displacement metrics.
#' @export
summary.ascfcd <- function(object, ...) {
  cat("==================================================\n")
  cat(" ASC-FCD: Multidimensional Functional Restructuring \n")
  cat("==================================================\n\n")

  cat(sprintf("Number of Contrasts: %d\n", length(object$rDelta_C)))
  cat(sprintf("Functional Dimensionality (k): %d (Explaining %.1f%% of variance)\n",
              object$k_retained, object$var_retained * 100))

  df_dist <- data.frame(
    Contrast = names(object$rDelta_C),
    Layer1_rDeltaC = round(object$rDelta_C, 2),
    Layer2_DeltaFDis = round(sapply(object$site_results, function(x) x$Delta_FDis), 4),
    Layer3_DeltaFRic = round(sapply(object$site_results, function(x) x$Delta_FRic), 4)
  )

  df_dist$Topology_Trend <- ifelse(df_dist$Layer3_DeltaFRic > 0, "Volume Expansion",
                                   ifelse(df_dist$Layer3_DeltaFRic < 0, "Volume Contraction",
                                          ifelse(df_dist$Layer2_DeltaFDis > 0.01, "Internal Expansion",
                                                 ifelse(df_dist$Layer2_DeltaFDis < -0.01, "Internal Contraction", "Stable"))))

  cat("\n--- Functional Shift Overview ---\n")
  print(df_dist, row.names = FALSE)
  cat("\n")

  if (!is.null(object$null_models)) {
    cat("--- Multi-Level Null Model Evaluation ---\n")
    cat("Struct: Incidence Filter | Quant: Demographic Filter | Identity: Trait Filter\n\n")
    print(object$null_models, row.names = FALSE)
    cat("\n")
  } else {
    cat("Null models not evaluated. Run asc_null() to test for structural, quantitative, and identity shifts.\n\n")
  }

  # TOP 5 DRIVERS (Leverage Preview)
  cat("--- Top Functional Drivers (Leverage Preview) ---\n")
  trans <- asc_transitions(object)
  df_drivers <- data.frame(Contrast = character(), Top_5_Drivers = character(), stringsAsFactors = FALSE)

  for(s in names(trans)) {
    top5 <- utils::head(trans[[s]]$species_leverage, 5)
    drivers_str <- paste(sprintf("%s (%+.2f)", top5$Species, top5$Leverage), collapse = ", ")
    df_drivers <- rbind(df_drivers, data.frame(Contrast = s, Top_5_Drivers = drivers_str))
  }
  print(df_drivers, row.names = FALSE)
  cat("\n")

  invisible(df_dist)
}


#' Summary for Pairwise Functional Spatial Divergence
#'
#' @description
#' Provides a comprehensive overview of the spatial network ASC-FCD analysis.
#'
#' @param object An object of class \code{ascfcd_pw}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the pairwise results data frame.
#' @export
summary.ascfcd_pw <- function(object, ...) {
  cat("==================================================\n")
  cat(" ASC-FCD: Pairwise Spatial Functional Network \n")
  cat("==================================================\n\n")

  n_com <- nrow(object$cwm_global)
  n_pairs <- nrow(object$pairwise_results)

  cat(sprintf("Communities Analyzed: %d | Spatial Contrasts: %d\n", n_com, n_pairs))
  cat(sprintf("Functional Dimensionality (k): %d (Explaining %.1f%% of variance)\n\n",
              object$k_retained, object$var_retained * 100))

  cat("--- Global Divergence Summary ---\n")
  rdelta <- object$pairwise_results$rDelta_C_pct
  fdis_shifts <- object$pairwise_results$Delta_FDis
  fric_shifts <- object$pairwise_results$Delta_FRic

  cat(sprintf(" Mean Position Shift (rDelta_C): %6.2f%%\n", mean(rdelta)))
  cat(sprintf(" Mean Dispersion Shift (Delta_FDis): %6.4f\n", mean(fdis_shifts)))
  cat(sprintf(" Mean Volume Shift (Delta_FRic): %6.4f\n\n", mean(fric_shifts)))

  if (!is.null(object$null_models)) {
    cat("--- Multi-Level Null Model Evaluation (Head) ---\n")
    cat("Struct: Incidence Filter | Quant: Demographic Filter | Identity: Trait Filter\n\n")
    print(utils::head(object$null_models, 10), row.names = FALSE)
    cat("\n")
  }

  # TOP 5 DRIVERS ESPACIALES (Leverage Preview)
  cat("--- Top Functional Drivers (Leverage Preview - Head) ---\n")
  trans <- asc_transitions(object)
  df_drivers <- data.frame(Contrast = character(), Top_5_Drivers = character(), stringsAsFactors = FALSE)

  for(s in names(trans)) {
    top5 <- utils::head(trans[[s]]$species_leverage, 5)
    drivers_str <- paste(sprintf("%s (%+.2f)", top5$Species, top5$Leverage), collapse = ", ")
    df_drivers <- rbind(df_drivers, data.frame(Contrast = s, Top_5_Drivers = drivers_str))
  }
  print(utils::head(df_drivers, 5), row.names = FALSE)
  cat("\n")

  invisible(object$pairwise_results)
}


#' Summary for Baseline Functional Entities
#'
#' @description
#' Provides an overview of the baseline functional topology of the evaluated entities.
#'
#' @param object An object of class \code{ascfcd_ent}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the entities results data frame.
#' @export
summary.ascfcd_base <- function(object, ...) {
  cat("==================================================\n")
  cat(" ASC-FCD: Baseline Functional Topology \n")
  cat("==================================================\n\n")

  cat(sprintf("Entities Analyzed: %d\n", nrow(object$entities_results)))
  cat(sprintf("Functional Dimensionality (k): %d (Explaining %.1f%% of variance)\n\n",
              object$k_retained, object$var_retained * 100))

  cat("--- Absolute Entity Metrics ---\n")

  # Redondeo para limpieza visual en consola
  print_df <- object$entities_results
  num_cols <- sapply(print_df, is.numeric)
  print_df[, num_cols] <- round(print_df[, num_cols], 4)

  print(print_df, row.names = FALSE)
  cat("\n")

  invisible(object$entities_results)
}

#' Summary for Functional Entities Classification
#'
#' @description
#' Provides an overview of the species clustering into functional entities.
#'
#' @param object An object of class \code{ascfcd_entities}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the classification data frame.
#' @export
summary.ascfcd_entities <- function(object, ...) {
  cat("==================================================\n")
  cat(" ASC-FCD: Functional Entities Classification \n")
  cat("==================================================\n\n")

  df <- object$entity_classification
  n_spp <- nrow(df)
  n_ent <- length(unique(df$Entity_ID))

  cat(sprintf("Total Species: %d\n", n_spp))
  cat(sprintf("Functional Entities Identified: %d\n", n_ent))
  cat(sprintf("Distance Metric: %s | Clustering Method: %s\n\n",
              tools::toTitleCase(object$metric), object$method))

  cat("--- Species Distribution per Entity ---\n")
  dist_table <- table(df$Entity_ID)
  print(dist_table)
  cat("\n")

  cat("--- Classification Preview (Head) ---\n")
  print(utils::head(df, 5), row.names = FALSE)
  cat("\n")

  invisible(df)
}
