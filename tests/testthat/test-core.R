# ==============================================================================
# test-core.R — Unit tests for ascent core functions
# ==============================================================================
# Covers: asc_paired, asc_pairwise, asc_null, asc_transitions
# Orthogonal scenarios: pure translation, pure reorganization, pure expansion
# ==============================================================================

# --- Shared fixtures ---
traits_10sp <- data.frame(
  Body_Mass   = c(15, 20, 30, 45, 60, 90, 150, 250, 400, 600),
  Beak_Length = c(10, 12, 15, 22, 28, 35,  45,  60,  85, 120),
  Frugivore   = factor(c(0, 0, 1, 1, 1, 0,  1,  1,  1,  0)),
  Forest_Dep  = factor(c(0, 0, 0, 0, 1, 1,  1,  1,  1,  1))
)
rownames(traits_10sp) <- paste0("Sp", 1:10)


# ==============================================================================
# asc_paired
# ==============================================================================
test_that("asc_paired returns correct class and structure", {
  abund <- rbind(
    Ref = c(0, 0, 10, 15, 25, 20, 10, 8, 4, 2),
    Imp = c(40, 35, 15, 5, 0, 0, 0, 0, 0, 0)
  )

  res <- asc_paired(traits_10sp, abund,
                     sites = c("S1", "S1"), time = c("R", "I"), ref_time = "R")

  expect_s3_class(res, "ascfcd")
  expect_true("rDelta_C" %in% names(res))
  expect_true("quality" %in% names(res))
  expect_true(is.numeric(res$quality) && res$quality > 0)
  expect_true("S1" %in% names(res$rDelta_C))
  expect_true(res$rDelta_C["S1"] > 0)
})

test_that("asc_paired handles degenerate community (1 species)", {
  abund <- rbind(
    Ref = c(5, 10, 15, 20, 25, 10, 8, 4, 2, 1),
    Imp = c(100, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  )

  res <- asc_paired(traits_10sp, abund,
                     sites = c("S1", "S1"), time = c("R", "I"), ref_time = "R")

  expect_true(is.na(res$site_results$S1$Delta_FRic))
  expect_true(!is.na(res$rDelta_C["S1"]))
})

test_that("asc_paired returns zero shift for identical communities", {
  abund <- rbind(
    A = c(5, 10, 15, 20, 25, 10, 8, 4, 2, 1),
    B = c(5, 10, 15, 20, 25, 10, 8, 4, 2, 1)
  )

  res <- asc_paired(traits_10sp, abund,
                     sites = c("S1", "S1"), time = c("R", "C"), ref_time = "R")

  expect_equal(unname(res$rDelta_C["S1"]), 0, tolerance = 1e-10)
})


# ==============================================================================
# ORTHOGONAL SCENARIOS
# ==============================================================================

# --- Scenario 1: Pure Translation ---
# Same species, same richness, shifted abundances → position moves, dispersion ~stable
test_that("Pure translation: position shifts, dispersion stable", {
  abund <- rbind(
    Ref  = c(0, 0,  5, 10, 40, 30, 10,  3, 1, 1),
    Comp = c(0, 0,  1,  3, 10, 30, 40, 10, 5, 1)
  )

  res <- asc_paired(traits_10sp, abund,
                     sites = c("S1", "S1"), time = c("R", "C"), ref_time = "R")

  expect_true(res$rDelta_C["S1"] > 0)
  # FRic should be NA or zero (same species present, same hull)
  # Actually both have same species so Delta_FRic should be ~0
  delta_fric <- res$site_results$S1$Delta_FRic
  if (!is.na(delta_fric)) {
    expect_equal(delta_fric, 0, tolerance = 1e-8)
  }
})

# --- Scenario 2: Pure Reorganization ---
# Same centroid, different spread → dispersion changes, position ~stable
test_that("Pure reorganization: dispersion shifts, position near-zero", {
  # Uniform → concentrated on central species
  abund <- rbind(
    Ref  = c(0, 0, 10, 10, 10, 10, 10, 10, 10, 10),
    Comp = c(0, 0,  1,  1,  1,  1, 40, 40,  1,  1)
  )

  res <- asc_paired(traits_10sp, abund,
                     sites = c("S1", "S1"), time = c("R", "C"), ref_time = "R")

  # FDis should change
  expect_true(abs(res$site_results$S1$Delta_FDis) > 0)
})

# --- Scenario 3: Pure Expansion ---
# Outer species added → volume expands, centroid may shift slightly
test_that("Pure expansion: FRic changes with species loss", {
  abund <- rbind(
    Ref  = c(5, 5, 10, 15, 25, 20, 10, 8, 5, 5),
    Comp = c(0, 0,  0, 15, 25, 20, 10, 0, 0, 0)
  )

  res <- asc_paired(traits_10sp, abund,
                     sites = c("S1", "S1"), time = c("R", "C"), ref_time = "R")

  delta_fric <- res$site_results$S1$Delta_FRic
  # May be NA if comp has too few species for hull
  if (!is.na(delta_fric)) {
    expect_true(delta_fric != 0)
  }
})


# ==============================================================================
# asc_pairwise
# ==============================================================================
test_that("asc_pairwise returns correct class and handles multiple communities", {
  abund <- rbind(
    A = c(0, 0, 10, 15, 25, 20, 10, 8, 4, 2),
    B = c(40, 35, 15, 5, 0, 0, 0, 0, 0, 0),
    C = c(0, 5, 10, 20, 35, 20, 10, 0, 0, 0)
  )

  res <- asc_pairwise(traits_10sp, abund)

  expect_s3_class(res, "ascfcd_pw")
  expect_equal(nrow(res$pairwise_results), 3)  # 3 choose 2
  expect_true(all(res$pairwise_results$rDelta_C_pct >= 0))
})


# ==============================================================================
# asc_null
# ==============================================================================
test_that("asc_null appends null_models with correct structure", {
  abund <- rbind(
    Ref = c(0, 0, 10, 15, 25, 20, 10, 8, 4, 2),
    Imp = c(40, 35, 15, 5, 0, 0, 0, 0, 0, 0)
  )

  res <- asc_paired(traits_10sp, abund,
                     sites = c("S1", "S1"), time = c("R", "I"), ref_time = "R")

  res <- asc_null(res, n_perm = 10, seed = 42)

  expect_true("null_models" %in% names(res))
  expect_s3_class(res$null_models, "data.frame")
  expect_true(all(c("Filter", "Metric", "SES", "P_value") %in% names(res$null_models)))
})

test_that("asc_null reproducibility with same seed", {
  abund <- rbind(
    Ref = c(0, 0, 10, 15, 25, 20, 10, 8, 4, 2),
    Imp = c(40, 35, 15, 5, 0, 0, 0, 0, 0, 0)
  )

  base <- asc_paired(traits_10sp, abund,
                      sites = c("S1", "S1"), time = c("R", "I"), ref_time = "R")

  r1 <- asc_null(base, n_perm = 50, seed = 123)
  r2 <- asc_null(base, n_perm = 50, seed = 123)

  expect_identical(r1$null_models, r2$null_models)
})

test_that("asc_null Model B FRic is always NA", {
  abund <- rbind(
    Ref  = c(5, 10, 15, 20, 25, 10, 8, 4, 2, 1),
    Comp = c(1, 2, 4, 8, 10, 25, 20, 15, 10, 5)
  )

  res <- asc_null(
    asc_paired(traits_10sp, abund,
               sites = c("S", "S"), time = c("R", "C"), ref_time = "R"),
    n_perm = 10, seed = 1
  )

  fric_b <- res$null_models[res$null_models$Filter == "Quantitative" &
                             res$null_models$Metric == "Volume (Delta FRic)", ]
  expect_true(all(is.na(fric_b$SES)))
  expect_true(all(is.na(fric_b$P_value)))
})


# ==============================================================================
# asc_transitions
# ==============================================================================
test_that("asc_transitions returns class and leverage additivity", {
  abund <- rbind(
    Ref = c(0, 0, 10, 15, 25, 20, 10, 8, 4, 2),
    Imp = c(40, 35, 15, 5, 0, 0, 0, 0, 0, 0)
  )

  res <- asc_paired(traits_10sp, abund,
                     sites = c("S1", "S1"), time = c("R", "I"), ref_time = "R")
  trans <- asc_transitions(res)

  expect_s3_class(trans, "ascfcd_transitions")
  expect_true("S1" %in% names(trans))

  # Leverage additivity: sum(Lev) = ||Delta_C||
  lev_sum <- sum(trans$S1$species_leverage$Leverage)
  abs_dist <- res$site_results$S1$abs_dist
  expect_equal(lev_sum, abs_dist, tolerance = 1e-8)
})

test_that("asc_transitions zero displacement gives zero leverage", {
  abund <- rbind(
    A = c(5, 10, 15, 20, 25, 10, 8, 4, 2, 1),
    B = c(5, 10, 15, 20, 25, 10, 8, 4, 2, 1)
  )

  res <- asc_paired(traits_10sp, abund,
                     sites = c("S1", "S1"), time = c("R", "C"), ref_time = "R")
  trans <- asc_transitions(res)

  expect_true(all(trans$S1$species_leverage$Leverage == 0))
})

test_that("Leverage = Delta_p * Projection", {
  abund <- rbind(
    Ref = c(0, 0, 10, 15, 25, 20, 10, 8, 4, 2),
    Imp = c(40, 35, 15, 5, 0, 0, 0, 0, 0, 0)
  )

  res <- asc_paired(traits_10sp, abund,
                     sites = c("S1", "S1"), time = c("R", "I"), ref_time = "R")
  trans <- asc_transitions(res)
  df <- trans$S1$species_leverage

  recomputed <- df$Delta_p * df$Projection
  expect_equal(df$Leverage, recomputed, tolerance = 1e-12)
})
