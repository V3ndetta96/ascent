# Evitar "Notes" en CRAN por variables globales de ggplot2
utils::globalVariables(c("Dim1", "Dim2", "Axis", "ASC", "Critical", "Type"))

#' Plot ASC-FCD Paired Contrast
#'
#' @description
#' Generates a publication-ready, 2-panel figure for a specific paired contrast.
#' The left panel shows the multidimensional functional space and the centroid trajectory.
#' The right panel displays the Axis-Specific Contribution (ASC) barplot, highlighting
#' the critical dimension of environmental change.
#'
#' @param x An object of class \code{ascfcd}.
#' @param contrast Character string specifying which site/community contrast to plot.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A \code{patchwork} object containing the combined ggplot2 figures.
#' @export
#'
#' @import ggplot2
#' @importFrom patchwork wrap_plots
plot.ascfcd <- function(x, contrast, ...) {

  if(missing(contrast) || !(contrast %in% names(x$rDelta_C))) {
    stop(sprintf("Please provide a valid contrast name. Available contrasts are: %s",
                 paste(names(x$rDelta_C), collapse = ", ")))
  }

  # ============================================================================
  # DATA PREPARATION
  # ============================================================================
  eje_critico <- x$critical_axes[contrast]

  # Si el eje crítico es el 1, graficamos 1 vs 2. Si es otro, graficamos 1 vs Crítico.
  eje_y <- ifelse(eje_critico == 1, 2, eje_critico)

  F_mat <- as.data.frame(x$F_matrix)
  colnames(F_mat) <- paste0("Dim", 1:ncol(F_mat))

  # Coordenadas de las especies (Fondo)
  df_spp <- data.frame(
    Dim1 = F_mat[, 1],
    Dim2 = F_mat[, eje_y]
  )

  # Coordenadas de los centroides (Trayectoria)
  c_ref <- x$CWM_ref[contrast, ]
  c_comp <- x$CWM_comp[contrast, ]

  df_traj <- data.frame(
    Dim1 = c(c_ref[1], c_comp[1]),
    Dim2 = c(c_ref[eje_y], c_comp[eje_y]),
    Type = c("Reference", "Comparison")
  )

  # Datos para el gráfico de barras ASC
  asc_vals <- x$ASC_j_list[[contrast]]
  df_asc <- data.frame(
    Axis = factor(1:length(asc_vals)),
    ASC = asc_vals,
    Critical = ifelse(1:length(asc_vals) == eje_critico, "Yes", "No")
  )

  # ============================================================================
  # PANEL A: FUNCTIONAL SPACE (PCoA)
  # ============================================================================
  var_x <- round(x$Var_j[1] * 100, 1)
  var_y <- round(x$Var_j[eje_y] * 100, 1)

  p1 <- ggplot() +
    # Nube de especies de fondo
    geom_point(data = df_spp, aes(x = Dim1, y = Dim2),
               color = "grey80", size = 1.5, alpha = 0.6) +
    # Punto de Referencia
    geom_point(data = subset(df_traj, Type == "Reference"),
               aes(x = Dim1, y = Dim2), color = "#3B9AB2", size = 4) +
    # Vector de trayectoria
    geom_path(data = df_traj, aes(x = Dim1, y = Dim2),
              color = "black", linewidth = 1,
              arrow = arrow(type = "closed", length = unit(0.15, "inches"))) +
    # Estética
    theme_bw(base_size = 14) +
    labs(
      title = sprintf("Functional Trajectory: %s", contrast),
      subtitle = sprintf("Relative Displacement (rDelta C): %.1f%%", x$rDelta_C[contrast]),
      x = sprintf("PCoA 1 (%.1f%%)", var_x),
      y = sprintf("PCoA %d (%.1f%%)", eje_y, var_y)
    ) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )

  # ============================================================================
  # PANEL B: ASC DECOMPOSITION
  # ============================================================================
  p2 <- ggplot(df_asc, aes(x = Axis, y = ASC, fill = Critical)) +
    geom_col(color = "black", width = 0.7) +
    scale_fill_manual(values = c("No" = "grey70", "Yes" = "#F21A00")) +
    theme_classic(base_size = 14) +
    labs(
      title = "Dimensional Drivers",
      subtitle = "Axis-Specific Contribution (ASC)",
      x = "Functional Dimension",
      y = "Contribution (%)"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )

  # ============================================================================
  # COMBINE AND RETURN
  # ============================================================================
  combined_plot <- patchwork::wrap_plots(p1, p2, widths = c(2, 1))
  return(combined_plot)
}
