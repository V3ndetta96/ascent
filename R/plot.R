#' Plot Multidimensional Functional Restructuring
#'
#' @description
#' Visualizes the functional trajectory (Centroid shift) and volume shift (Convex Hull)
#' of a specific site/contrast in the PCoA space, alongside the species leverage.
#'
#' @param x An object of class \code{ascfcd}.
#' @param contrast Character. The exact name of the site/contrast to plot.
#' @param type Character. What to plot? \code{"pcoa"} (Trajectory & Hull), \code{"leverage"} (Drivers), or \code{"both"}. Default is \code{"both"}.
#' @param n_sp Integer. Number of top species to display in the leverage plot. Default is 10.
#' @param ... Additional graphical arguments.
#'
#' @import ggplot2
#' @return A ggplot object (or patchwork object if \code{type = "both"}).
#' @export
plot.ascfcd <- function(x, contrast, type = c("both", "pcoa", "leverage"), n_sp = 10, ...) {

  type <- match.arg(type)
  if(!(contrast %in% names(x$rDelta_C))) {
    stop(sprintf("Contrast '%s' not found. Available: %s", contrast, paste(names(x$rDelta_C), collapse = ", ")))
  }

  # ============================================================================
  # 1. PCOA TRAJECTORY & HULL PLOT
  # ============================================================================
  F_mat <- x$F_space
  c_ref <- x$cwm_ref[[contrast]]
  c_comp <- x$cwm_comp[[contrast]]

  # Extracción de abundancias para este contraste
  p_ref <- as.numeric(x$p_ref[[contrast]])
  p_comp <- as.numeric(x$p_comp[[contrast]])
  sp_names <- rownames(F_mat)

  # Coordenadas de especies (Fondo)
  df_spp <- as.data.frame(F_mat[, 1:2, drop = FALSE])
  colnames(df_spp) <- c("PC1", "PC2")

  # Coordenadas de Centroides
  df_cent <- data.frame(
    PC1 = c(c_ref[1], c_comp[1]),
    PC2 = c(c_ref[2], c_comp[2]),
    State = factor(c("Reference", "Comparison"), levels = c("Reference", "Comparison"))
  )

  # Cálculo de polígonos Convex Hull (Capa 3)
  calc_hull_2d <- function(abund_vec, state_name) {
    idx <- abund_vec > 0
    if(sum(idx) < 3) return(data.frame(PC1=numeric(0), PC2=numeric(0), State=character(0)))
    coords <- df_spp[idx, ]
    hull_idx <- grDevices::chull(coords$PC1, coords$PC2)
    poly_df <- coords[hull_idx, ]
    poly_df$State <- state_name
    return(poly_df)
  }

  df_hulls <- rbind(
    calc_hull_2d(p_ref, "Reference"),
    calc_hull_2d(p_comp, "Comparison")
  )
  if(nrow(df_hulls) > 0) {
    df_hulls$State <- factor(df_hulls$State, levels = c("Reference", "Comparison"))
  }

  # Textos de varianza y diagnóstico
  var_p1 <- x$axis_var[1] * 100
  var_p2 <- x$axis_var[2] * 100

  res <- x$site_results[[contrast]]
  trend <- ifelse(res$Delta_FRic > 0, "Vol Expansion", ifelse(res$Delta_FRic < 0, "Vol Contraction",
                                                              ifelse(res$Delta_FDis > 0.01, "Internal Expansion", ifelse(res$Delta_FDis < -0.01, "Internal Contraction", "Stable"))))

  # Subtítulo con salto de línea (\n) para evitar superposiciones
  sub_text <- sprintf("Delta C: %.1f%%\nDelta FDis: %.2f | Trend: %s",
                      x$rDelta_C[contrast], res$Delta_FDis, trend)

  p_pca <- ggplot() +
    geom_point(data = df_spp, aes(x = PC1, y = PC2), color = "grey85", size = 1.5, alpha = 0.8)

  # Añadir polígonos si existen
  if(nrow(df_hulls) > 0) {
    p_pca <- p_pca + geom_polygon(data = df_hulls, aes(x = PC1, y = PC2, fill = State, color = State), alpha = 0.15, linewidth = 0.5)
  }

  p_pca <- p_pca +
    geom_segment(data = df_cent, aes(x = PC1[1], y = PC2[1], xend = PC1[2], yend = PC2[2]),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"), color = "black", linewidth = 0.8) +
    geom_point(data = df_cent, aes(x = PC1, y = PC2, color = State), size = 4) +
    scale_color_manual(values = c("Reference" = "#4575b4", "Comparison" = "#d73027")) +
    scale_fill_manual(values = c("Reference" = "#4575b4", "Comparison" = "#d73027")) +
    theme_minimal(base_size = 13) +
    labs(title = paste("Topology:", contrast), subtitle = sub_text,
         x = sprintf("PCoA 1 (%.1f%%)", var_p1), y = sprintf("PCoA 2 (%.1f%%)", var_p2)) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.subtitle = element_text(size = 10, lineheight = 1.2, color = "grey30"))

  if (type == "pcoa") return(p_pca)

  # ============================================================================
  # 2. FUNCTIONAL LEVERAGE PLOT
  # ============================================================================
  trans <- asc_transitions(x)
  df_lev <- trans[[contrast]]$species_leverage
  df_lev$Abs_Lev <- abs(df_lev$Leverage)
  df_lev <- df_lev[order(df_lev$Abs_Lev, decreasing = TRUE), ]
  df_lev <- utils::head(df_lev, n_sp)
  df_lev <- df_lev[order(df_lev$Leverage), ]
  df_lev$Species <- factor(df_lev$Species, levels = df_lev$Species)

  p_lev <- ggplot(df_lev, aes(x = Leverage, y = Species)) +
    geom_segment(aes(x = 0, xend = Leverage, y = Species, yend = Species), color = "grey60", linewidth = 1) +
    geom_point(aes(color = Leverage > 0), size = 4) +
    scale_color_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4")) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.8) +
    theme_minimal(base_size = 13) +
    labs(title = "Functional Leverage", subtitle = "Top Drivers of Centroid Shift",
         x = "Leverage Score", y = "") +
    theme(legend.position = "none", panel.grid.minor = element_blank())

  if (type == "leverage") return(p_lev)

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("Package 'patchwork' is required to plot type = 'both'. Returning PCoA only.")
    return(p_pca)
  }
  return(p_pca | p_lev)
}

#' Plot Pairwise Functional Spatial Divergence
#'
#' @description
#' Visualizes spatial functional restructuring (Centroid & Hull) between two communities.
#'
#' @param x An object of class \code{ascfcd_pw}.
#' @param contrast Character. The contrast to plot (e.g., "SiteA_vs_SiteB").
#' @param type Character. What to plot? \code{"pcoa"}, \code{"leverage"}, or \code{"both"}.
#' @param n_sp Integer. Number of top species to display. Default is 10.
#' @param ... Additional graphical arguments.
#'
#' @import ggplot2
#' @return A ggplot object.
#' @export
plot.ascfcd_pw <- function(x, contrast, type = c("both", "pcoa", "leverage"), n_sp = 10, ...) {

  type <- match.arg(type)
  if(!(contrast %in% names(x$directional_vectors))) {
    stop(sprintf("Contrast '%s' not found. Available: %s", contrast, paste(utils::head(names(x$directional_vectors)), collapse = ", ")))
  }

  F_mat <- x$F_space
  coms <- unlist(strsplit(contrast, "_vs_"))
  com_A <- coms[1]; com_B <- coms[2]
  c_A <- x$cwm_global[com_A, ]; c_B <- x$cwm_global[com_B, ]

  df_spp <- as.data.frame(F_mat[, 1:2, drop = FALSE])
  colnames(df_spp) <- c("PC1", "PC2")

  df_cent <- data.frame(
    PC1 = c(c_A[1], c_B[1]), PC2 = c(c_A[2], c_B[2]),
    Community = factor(c(com_A, com_B), levels = c(com_A, com_B))
  )

  calc_hull_2d <- function(abund_vec, state_name) {
    idx <- abund_vec > 0
    if(sum(idx) < 3) return(data.frame(PC1=numeric(0), PC2=numeric(0), Community=character(0)))
    coords <- df_spp[idx, ]
    hull_idx <- grDevices::chull(coords$PC1, coords$PC2)
    poly_df <- coords[hull_idx, ]
    poly_df$Community <- state_name
    return(poly_df)
  }

  df_hulls <- rbind(
    calc_hull_2d(as.numeric(x$p_ref[[contrast]]), com_A),
    calc_hull_2d(as.numeric(x$p_comp[[contrast]]), com_B)
  )
  if(nrow(df_hulls) > 0) df_hulls$Community <- factor(df_hulls$Community, levels = c(com_A, com_B))

  res <- x$pairwise_results[x$pairwise_results$Community_A == com_A & x$pairwise_results$Community_B == com_B, ]
  trend <- ifelse(res$Delta_FRic > 0, "Vol Expansion", ifelse(res$Delta_FRic < 0, "Vol Contraction",
                                                              ifelse(res$Delta_FDis > 0.01, "Internal Expansion", ifelse(res$Delta_FDis < -0.01, "Internal Contraction", "Stable"))))

  # Subtítulo con salto de línea (\n) para evitar superposiciones
  sub_text <- sprintf("rDelta C: %.1f%%\nDelta FDis: %.2f | Trend: %s", res$rDelta_C_pct, res$Delta_FDis, trend)

  p_pca <- ggplot() +
    geom_point(data = df_spp, aes(x = PC1, y = PC2), color = "grey85", size = 1.5, alpha = 0.8)

  if(nrow(df_hulls) > 0) {
    p_pca <- p_pca + geom_polygon(data = df_hulls, aes(x = PC1, y = PC2, fill = Community, color = Community), alpha = 0.15, linewidth = 0.5)
  }

  p_pca <- p_pca +
    geom_segment(data = df_cent, aes(x = PC1[1], y = PC2[1], xend = PC1[2], yend = PC2[2]),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"), color = "black", linewidth = 0.8) +
    geom_point(data = df_cent, aes(x = PC1, y = PC2, color = Community), size = 4) +
    scale_color_manual(values = stats::setNames(c("#4575b4", "#d73027"), c(com_A, com_B))) +
    scale_fill_manual(values = stats::setNames(c("#4575b4", "#d73027"), c(com_A, com_B))) +
    theme_minimal(base_size = 13) +
    labs(title = paste("Spatial Topology:", contrast), subtitle = sub_text,
         x = sprintf("PCoA 1 (%.1f%%)", x$axis_var[1]*100), y = sprintf("PCoA 2 (%.1f%%)", x$axis_var[2]*100)) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.subtitle = element_text(size = 10, lineheight = 1.2, color = "grey30"))

  if (type == "pcoa") return(p_pca)

  # Leverage Espacial
  delta_p <- as.numeric(x$p_comp[[contrast]] - x$p_ref[[contrast]])
  v_dir <- x$directional_vectors[[contrast]]
  mag_v <- sqrt(sum(v_dir^2))

  if(mag_v == 0) {
    df_lev <- data.frame(Species = rownames(F_mat), Leverage = 0)
  } else {
    u_dir <- v_dir / mag_v
    projections <- numeric(nrow(F_mat))
    for (i in seq_len(nrow(F_mat))) {
      projections[i] <- sum((F_mat[i, ] - as.numeric(c_A)) * u_dir)
    }
    df_lev <- data.frame(Species = rownames(F_mat), Leverage = delta_p * projections)
  }

  df_lev$Abs_Lev <- abs(df_lev$Leverage)
  df_lev <- df_lev[order(df_lev$Abs_Lev, decreasing = TRUE), ]
  df_lev <- utils::head(df_lev, n_sp)
  df_lev <- df_lev[order(df_lev$Leverage), ]
  df_lev$Species <- factor(df_lev$Species, levels = df_lev$Species)

  p_lev <- ggplot(df_lev, aes(x = Leverage, y = Species)) +
    geom_segment(aes(x = 0, xend = Leverage, y = Species, yend = Species), color = "grey60", linewidth = 1) +
    geom_point(aes(color = Leverage > 0), size = 4) +
    scale_color_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4")) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.8) +
    theme_minimal(base_size = 13) +
    labs(title = "Spatial Leverage", subtitle = "Drivers of divergence", x = "Leverage Score", y = "") +
    theme(legend.position = "none", panel.grid.minor = element_blank())

  if (type == "leverage") return(p_lev)
  if (!requireNamespace("patchwork", quietly = TRUE)) return(p_pca)
  return(p_pca | p_lev)
}
