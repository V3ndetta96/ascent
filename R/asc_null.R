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
#' @note
#' **Model A (Structural):** After curveball permutation, all species present
#' receive uniform relative abundance (1/S_local). This means the null distribution
#' for FDis conflates the effect of taxonomic identity with the assumption of
#' equitability. The SES tests whether the observed shift is extreme given
#' random species composition, *not* given random composition with the observed SAD.
#'
#' **Model B (Quantitative):** Because incidence is held fixed, the convex hull
#' is invariant across permutations. Delta FRic under Model B is reported as
#' \code{NA} (not applicable), not zero.
#'
#' **Model C (Identity):** Statistical power is limited when the regional species
#' pool is small (< 15 species). With few species, the number of unique trait
#' permutations is small, reducing the resolution of the null distribution.
#' Consider increasing \code{n_perm} and interpreting marginal p-values
#' (0.05 < p < 0.10) with caution.
#'
#' **Scope:** When multiple contrasts are evaluated simultaneously, the curveball
#' operates on the full stacked matrix, assuming a shared regional species pool.
#' For biogeographically independent sites, run \code{asc_null()} on each
#' contrast separately.
#'
#' @importFrom vegan nullmodel
#' @importFrom stats sd simulate
#' @return The original object with an appended \code{null_models} data frame.
#'
#' @examples
#' \donttest{
#' traits <- data.frame(
#'   Mass = c(15, 30, 60, 150, 400),
#'   Beak = c(10, 15, 28,  45,  85),
#'   Diet = factor(c(0, 1, 1, 1, 0))
#' )
#' rownames(traits) <- paste0("Sp", 1:5)
#'
#' abund <- rbind(
#'   Reference = c(0, 10, 25, 20, 5),
#'   Impacted  = c(30, 15,  5,  0, 0)
#' )
#'
#' res <- asc_paired(
#'   traits, abund,
#'   sites = c("S1", "S1"),
#'   time = c("Reference", "Impacted"),
#'   ref_time = "Reference"
#' )
#' res <- asc_null(res, n_perm = 99, seed = 42)
#' res$null_models
#' }
#'
#' @export
asc_null <- function(x, n_perm = 999, seed = NULL) {

  if(!inherits(x, c("ascfcd", "ascfcd_pw"))) stop("Input object must be of class 'ascfcd' or 'ascfcd_pw'.")

  # Safe RNG handling: restore .Random.seed on exit
  if(!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- .GlobalEnv$.Random.seed
      on.exit(.GlobalEnv$.Random.seed <- old_seed, add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }

  F_mat <- x$F_space
  entities <- names(x$directional_vectors)

  # Build M matrix with explicit entity-to-row mapping
  M_list <- list()
  entity_row_ref <- integer(length(entities))
  entity_row_comp <- integer(length(entities))
  names(entity_row_ref) <- names(entity_row_comp) <- entities

  for(s in entities) {
    M_list[[paste0(s, "_ref")]] <- x$p_ref[[s]]
    M_list[[paste0(s, "_comp")]] <- x$p_comp[[s]]
  }
  M <- do.call(rbind, M_list)

  # Map each entity to its ref/comp row indices in M (explicit, not positional)
  for (j in seq_along(entities)) {
    entity_row_ref[entities[j]]  <- 2L * j - 1L
    entity_row_comp[entities[j]] <- 2L * j
  }

  # Validation: row count must equal 2 * number of entities
  stopifnot(nrow(M) == 2L * length(entities))

  M_bin <- ifelse(M > 0, 1, 0)
  nm_struct <- vegan::nullmodel(M_bin, method = "curveball")
  sims_struct <- stats::simulate(nm_struct, nsim = n_perm)

  shuffle_nonzero <- function(vec) {
    idx <- vec > 0
    if(sum(idx) > 1) vec[idx] <- sample(vec[idx])
    return(vec)
  }

  # Containers for Model A (Structural) and Model C (Identity)
  null_pos_A <- matrix(NA, n_perm, length(entities)); colnames(null_pos_A) <- entities
  null_pos_B <- matrix(NA, n_perm, length(entities)); colnames(null_pos_B) <- entities
  null_pos_C <- matrix(NA, n_perm, length(entities)); colnames(null_pos_C) <- entities

  null_fdis_A <- matrix(NA, n_perm, length(entities)); colnames(null_fdis_A) <- entities
  null_fdis_B <- matrix(NA, n_perm, length(entities)); colnames(null_fdis_B) <- entities
  null_fdis_C <- matrix(NA, n_perm, length(entities)); colnames(null_fdis_C) <- entities

  null_fric_A <- matrix(NA, n_perm, length(entities)); colnames(null_fric_A) <- entities
  # null_fric_B intentionally omitted: incidence is fixed under Model B,
  # so the convex hull is invariant. FRic has no null distribution.
  null_fric_C <- matrix(NA, n_perm, length(entities)); colnames(null_fric_C) <- entities

  for (i in 1:n_perm) {
    # Simulated matrices (A and B)
    mat_A_bin <- sims_struct[, , i]
    mat_A_rel <- sweep(mat_A_bin, 1, rowSums(mat_A_bin), "/"); mat_A_rel[is.na(mat_A_rel)] <- 0
    # M already contains relative abundances; shuffle preserves row sums.
    # Re-normalization is defensive (guards against floating-point drift).
    mat_B <- t(apply(M, 1, shuffle_nonzero))
    mat_B_rel <- mat_B  # already relative; no sweep needed

    # Model C: Shuffled Trait Space
    F_mat_shuf <- F_mat[sample(nrow(F_mat)), , drop = FALSE]

    row_idx <- 1
    for(s in entities) {
      r_idx <- entity_row_ref[s]
      c_idx <- entity_row_comp[s]

      # --- MODEL A: STRUCTURAL (Curveball) ---
      c_r_A <- mat_A_rel[r_idx, , drop=FALSE] %*% F_mat
      c_c_A <- mat_A_rel[c_idx, , drop=FALSE] %*% F_mat
      null_pos_A[i, s] <- sqrt(sum((c_c_A - c_r_A)^2))
      topo_r_A <- .calc_core_topology(mat_A_rel[r_idx, ], F_mat)
      topo_c_A <- .calc_core_topology(mat_A_rel[c_idx, ], F_mat)
      null_fdis_A[i, s] <- topo_c_A$FDis - topo_r_A$FDis
      null_fric_A[i, s] <- topo_c_A$FRic - topo_r_A$FRic

      # --- MODEL B: QUANTITATIVE (SAD Reshuffle) ---
      c_r_B <- mat_B_rel[r_idx, , drop=FALSE] %*% F_mat
      c_c_B <- mat_B_rel[c_idx, , drop=FALSE] %*% F_mat
      null_pos_B[i, s] <- sqrt(sum((c_c_B - c_r_B)^2))
      topo_r_B <- .calc_core_topology(mat_B_rel[r_idx, ], F_mat)
      topo_c_B <- .calc_core_topology(mat_B_rel[c_idx, ], F_mat)
      null_fdis_B[i, s] <- topo_c_B$FDis - topo_r_B$FDis
      # FRic under Model B: not evaluated (incidence fixed -> hull invariant)

      # --- MODEL C: IDENTITY (Trait Shuffle) ---
      c_r_C <- M[r_idx, , drop=FALSE] %*% F_mat_shuf
      c_c_C <- M[c_idx, , drop=FALSE] %*% F_mat_shuf
      null_pos_C[i, s] <- sqrt(sum((c_c_C - c_r_C)^2))
      topo_r_C <- .calc_core_topology(as.numeric(M[r_idx, ]), F_mat_shuf)
      topo_c_C <- .calc_core_topology(as.numeric(M[c_idx, ]), F_mat_shuf)
      null_fdis_C[i, s] <- topo_c_C$FDis - topo_r_C$FDis
      null_fric_C[i, s] <- topo_c_C$FRic - topo_r_C$FRic
    }
  }

  # Evaluation helper with robust sd guard
  eval_null <- function(obs, null_dist, two_tailed = FALSE) {
    if (all(is.na(null_dist))) return(list(SES = NA_real_, P = NA_real_))
    null_dist_clean <- null_dist[!is.na(null_dist)]
    if (length(null_dist_clean) < 2) return(list(SES = NA_real_, P = NA_real_))
    sd_val <- stats::sd(null_dist_clean)
    if (sd_val < 1e-12) return(list(SES = NaN, P = 1.000))
    ses <- (obs - mean(null_dist_clean)) / sd_val
    n_eff <- length(null_dist_clean)
    pval <- if(two_tailed) {
      (sum(abs(null_dist_clean) >= abs(obs)) + 1) / (n_eff + 1)
    } else {
      (sum(null_dist_clean >= obs) + 1) / (n_eff + 1)
    }
    return(list(SES = ses, P = pval))
  }

  df_res <- data.frame()
  for(s in entities) {
    obs_pos <- if(!is.null(x$site_results[[s]]$abs_dist)) x$site_results[[s]]$abs_dist else x$pairwise_results$Delta_C_abs[x$pairwise_results$Community_A == strsplit(s, "_vs_")[[1]][1] & x$pairwise_results$Community_B == strsplit(s, "_vs_")[[1]][2]]
    obs_fdis <- if(!is.null(x$site_results[[s]]$Delta_FDis)) x$site_results[[s]]$Delta_FDis else x$pairwise_results$Delta_FDis[x$pairwise_results$Community_A == strsplit(s, "_vs_")[[1]][1] & x$pairwise_results$Community_B == strsplit(s, "_vs_")[[1]][2]]
    obs_fric <- if(!is.null(x$site_results[[s]]$Delta_FRic)) x$site_results[[s]]$Delta_FRic else x$pairwise_results$Delta_FRic[x$pairwise_results$Community_A == strsplit(s, "_vs_")[[1]][1] & x$pairwise_results$Community_B == strsplit(s, "_vs_")[[1]][2]]

    e_pos_A <- eval_null(obs_pos, null_pos_A[, s], FALSE)
    e_pos_B <- eval_null(obs_pos, null_pos_B[, s], FALSE)
    e_pos_C <- eval_null(obs_pos, null_pos_C[, s], FALSE)

    e_fdis_A <- eval_null(obs_fdis, null_fdis_A[, s], TRUE)
    e_fdis_B <- eval_null(obs_fdis, null_fdis_B[, s], TRUE)
    e_fdis_C <- eval_null(obs_fdis, null_fdis_C[, s], TRUE)

    e_fric_A <- eval_null(obs_fric, null_fric_A[, s], TRUE)
    # FRic under Model B: bypass (not applicable)
    e_fric_B <- list(SES = NA_real_, P = NA_real_)
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
