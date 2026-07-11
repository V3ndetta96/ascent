# ==============================================================================
# test-assess.R — Unit tests for assess_functional_space
# ==============================================================================

test_that("assess_functional_space returns correct structure and values", {
  traits <- data.frame(
    Mass = c(15, 30, 60, 150, 400),
    Beak = c(10, 15, 28,  45,  85)
  )
  rownames(traits) <- paste0("Sp", 1:5)

  diag <- assess_functional_space(traits, dist_method = "euclidean")

  expect_s3_class(diag, "ascfcd_space")
  expect_equal(diag$n_species, 5)
  expect_equal(diag$n_traits, 2)
  expect_true(diag$quality > 0)
  expect_true(is.numeric(diag$eigenvalues))
})
