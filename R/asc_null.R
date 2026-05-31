#' Multi-Level Null Models for Multidimensional Functional Shifts
#'
#' @description
#' Evaluates the statistical significance of functional shifts using three null models:
#' Structural (Curveball), Quantitative (SAD Reshuffle), and Identity (Trait Shuffle).
#'
#' @param x An object of class \code{ascfcd} or \code{ascfcd_pw}.
#' @param n_perm Integer. Number of permutations. Default is 999.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @importFrom vegan nullmodel
#' @importFrom stats sd simulate
#' @return The original object with an appended \code{null_models} data frame.
#' @export
asc_null <- function(x, n_perm = 999, seed = NULL) {

  if(!inherits(x, c("ascfcd", "ascfcd_pw"))) stop("Input object must be of class 'ascfcd' or 'ascfcd_pw'.")
  if (!requireNamespace("vegan", quietly = TRUE)) stop("The 'vegan' package is required.")
  if (!requireNamespace("geometry", quietly = TRUE)) stop("The 'geometry' package is required.")

  if(!is.null(seed)) set.seed(seed)
  F_mat <- x$F_space
  entities <- names(x$directional_vectors)

  M_list <- list()
  for(s in entities) {
    M_list[[paste0(s, "_ref")]] <- x$p_ref[[s]]
    M_list[[paste0(s, "_comp")]] <- x$p_comp[[s]]
  }
  M <- do.call(rbind, M_list)

  M_bin <- ifelse(M > 0, 1, 0)
  nm_struct <- vegan::nullmodel(M_bin, method = "curveball")
  sims_struct <- stats::simulate(nm_struct, nsim = n_perm)

  shuffle_nonzero <- function(vec) {
    idx <- vec > 0
    if(sum(idx) > 1) vec[idx] <- sample(vec[idx])
    return(vec)
  }

  calc_topo <- function(p_vec, f_space) {
    idx <- p_vec > 0
    if(sum(idx) < 2) return(list(FDis=0, FRic=0))
    p_pres <- p_vec[idx] / sum(p_vec[idx])
    f_pres <- f_space[idx, , drop = FALSE]
    cent <- colSums(p_pres * f_pres)
    fdis <- sum(p_pres * sqrt(rowSums(sweep(f_pres, 2, cent, "-")^2)))

    fric <- 0
    if (sum(idx) > ncol(f_space)) {
      fric <- tryCatch({ geometry::convhulln(f_pres, options = "FA")$vol }, error = function(e) 0)
    }
    return(list(FDis = fdis, FRic = fric))
  }

  # Contenedores para las 3 Familias de Modelos Nulos
  null_pos_A <- matrix(NA, n_perm, length(entities)); colnames(null_pos_A) <- entities
  null_pos_B <- matrix(NA, n_perm, length(entities)); colnames(null_pos_B) <- entities
  null_pos_C <- matrix(NA, n_perm, length(entities)); colnames(null_pos_C) <- entities

  null_fdis_A <- matrix(NA, n_perm, length(entities)); colnames(null_fdis_A) <- entities
  null_fdis_B <- matrix(NA, n_perm, length(entities)); colnames(null_fdis_B) <- entities
  null_fdis_C <- matrix(NA, n_perm, length(entities)); colnames(null_fdis_C) <- entities

  null_fric_A <- matrix(NA, n_perm, length(entities)); colnames(null_fric_A) <- entities
  null_fric_B <- matrix(NA, n_perm, length(entities)); colnames(null_fric_B) <- entities
  null_fric_C <- matrix(NA, n_perm, length(entities)); colnames(null_fric_C) <- entities

  for (i in 1:n_perm) {
    # Matrices simuladas (A y B)
    mat_A_bin <- sims_struct[, , i]
    mat_A_rel <- sweep(mat_A_bin, 1, rowSums(mat_A_bin), "/"); mat_A_rel[is.na(mat_A_rel)] <- 0
    mat_B <- t(apply(M, 1, shuffle_nonzero))
    mat_B_rel <- sweep(mat_B, 1, rowSums(mat_B), "/")

    # Modelo C: Shuffled Trait Space (Identidades funcionales aleatorias)
    F_mat_shuf <- F_mat[sample(nrow(F_mat)), , drop = FALSE]

    row_idx <- 1
    for(s in entities) {

      # --- MODELO A: ESTRUCTURAL (Curveball) ---
      c_r_A <- mat_A_rel[row_idx, , drop=FALSE] %*% F_mat
      c_c_A <- mat_A_rel[row_idx+1, , drop=FALSE] %*% F_mat
      null_pos_A[i, s] <- sqrt(sum((c_c_A - c_r_A)^2))
      topo_r_A <- calc_topo(mat_A_rel[row_idx, ], F_mat)
      topo_c_A <- calc_topo(mat_A_rel[row_idx+1, ], F_mat)
      null_fdis_A[i, s] <- topo_c_A$FDis - topo_r_A$FDis
      null_fric_A[i, s] <- topo_c_A$FRic - topo_r_A$FRic

      # --- MODELO B: CUANTITATIVO (SAD Reshuffle) ---
      c_r_B <- mat_B_rel[row_idx, , drop=FALSE] %*% F_mat
      c_c_B <- mat_B_rel[row_idx+1, , drop=FALSE] %*% F_mat
      null_pos_B[i, s] <- sqrt(sum((c_c_B - c_r_B)^2))
      topo_r_B <- calc_topo(mat_B_rel[row_idx, ], F_mat)
      topo_c_B <- calc_topo(mat_B_rel[row_idx+1, ], F_mat)
      null_fdis_B[i, s] <- topo_c_B$FDis - topo_r_B$FDis
      null_fric_B[i, s] <- topo_c_B$FRic - topo_r_B$FRic

      # --- MODELO C: IDENTIDAD (Trait Shuffle / Riqueza Fija) ---
      # Usamos las abundancias empíricas reales, pero el espacio funcional permutado
      mat_r_real <- M[row_idx, , drop=FALSE] / sum(M[row_idx, ])
      mat_c_real <- M[row_idx+1, , drop=FALSE] / sum(M[row_idx+1, ])

      c_r_C <- mat_r_real %*% F_mat_shuf
      c_c_C <- mat_c_real %*% F_mat_shuf
      null_pos_C[i, s] <- sqrt(sum((c_c_C - c_r_C)^2))
      topo_r_C <- calc_topo(mat_r_real, F_mat_shuf)
      topo_c_C <- calc_topo(mat_c_real, F_mat_shuf)
      null_fdis_C[i, s] <- topo_c_C$FDis - topo_r_C$FDis
      null_fric_C[i, s] <- topo_c_C$FRic - topo_r_C$FRic

      row_idx <- row_idx + 2
    }
  }

  df_res <- data.frame()
  for(s in entities) {
    obs_pos <- if(!is.null(x$site_results[[s]]$abs_dist)) x$site_results[[s]]$abs_dist else x$pairwise_results$Delta_C_abs[x$pairwise_results$Community_A == strsplit(s, "_vs_")[[1]][1] & x$pairwise_results$Community_B == strsplit(s, "_vs_")[[1]][2]]
    obs_fdis <- if(!is.null(x$site_results[[s]]$Delta_FDis)) x$site_results[[s]]$Delta_FDis else x$pairwise_results$Delta_FDis[x$pairwise_results$Community_A == strsplit(s, "_vs_")[[1]][1] & x$pairwise_results$Community_B == strsplit(s, "_vs_")[[1]][2]]
    obs_fric <- if(!is.null(x$site_results[[s]]$Delta_FRic)) x$site_results[[s]]$Delta_FRic else x$pairwise_results$Delta_FRic[x$pairwise_results$Community_A == strsplit(s, "_vs_")[[1]][1] & x$pairwise_results$Community_B == strsplit(s, "_vs_")[[1]][2]]

    eval_null <- function(obs, null_dist, two_tailed = FALSE) {
      if(stats::sd(null_dist) == 0) return(list(SES = NaN, P = 1.000))
      ses <- (obs - mean(null_dist)) / stats::sd(null_dist)
      pval <- if(two_tailed) (sum(abs(null_dist) >= abs(obs)) + 1) / (n_perm + 1) else (sum(null_dist >= obs) + 1) / (n_perm + 1)
      return(list(SES = ses, P = pval))
    }

    e_pos_A <- eval_null(obs_pos, null_pos_A[, s], FALSE)
    e_pos_B <- eval_null(obs_pos, null_pos_B[, s], FALSE)
    e_pos_C <- eval_null(obs_pos, null_pos_C[, s], FALSE)

    e_fdis_A <- eval_null(obs_fdis, null_fdis_A[, s], TRUE)
    e_fdis_B <- eval_null(obs_fdis, null_fdis_B[, s], TRUE)
    e_fdis_C <- eval_null(obs_fdis, null_fdis_C[, s], TRUE)

    e_fric_A <- eval_null(obs_fric, null_fric_A[, s], TRUE)
    e_fric_B <- eval_null(obs_fric, null_fric_B[, s], TRUE)
    e_fric_C <- eval_null(obs_fric, null_fric_C[, s], TRUE)

    df_res <- rbind(df_res,
                    data.frame(Contrast = s, Filter = "Structural", Metric = "Position (Delta C)", Observed = round(obs_pos, 4), SES = round(e_pos_A$SES, 2), P_value = round(e_pos_A$P, 3)),
                    data.frame(Contrast = s, Filter = "Quantitative", Metric = "Position (Delta C)", Observed = round(obs_pos, 4), SES = round(e_pos_B$SES, 2), P_value = round(e_pos_B$P, 3)),
                    data.frame(Contrast = s, Filter = "Identity", Metric = "Position (Delta C)", Observed = round(obs_pos, 4), SES = round(e_pos_C$SES, 2), P_value = round(e_pos_C$P, 3)),

                    data.frame(Contrast = s, Filter = "Structural", Metric = "Dispersion (Delta FDis)", Observed = round(obs_fdis, 4), SES = round(e_fdis_A$SES, 2), P_value = round(e_fdis_A$P, 3)),
                    data.frame(Contrast = s, Filter = "Quantitative", Metric = "Dispersion (Delta FDis)", Observed = round(obs_fdis, 4), SES = round(e_fdis_B$SES, 2), P_value = round(e_fdis_B$P, 3)),
                    data.frame(Contrast = s, Filter = "Identity", Metric = "Dispersion (Delta FDis)", Observed = round(obs_fdis, 4), SES = round(e_fdis_C$SES, 2), P_value = round(e_fdis_C$P, 3)),

                    data.frame(Contrast = s, Filter = "Structural", Metric = "Volume (Delta FRic)", Observed = round(obs_fric, 4), SES = round(e_fric_A$SES, 2), P_value = round(e_fric_A$P, 3)),
                    data.frame(Contrast = s, Filter = "Quantitative", Metric = "Volume (Delta FRic)", Observed = round(obs_fric, 4), SES = round(e_fric_B$SES, 2), P_value = round(e_fric_B$P, 3)),
                    data.frame(Contrast = s, Filter = "Identity", Metric = "Volume (Delta FRic)", Observed = round(obs_fric, 4), SES = round(e_fric_C$SES, 2), P_value = round(e_fric_C$P, 3))
    )
  }
  x$null_models <- df_res
  return(x)
}
