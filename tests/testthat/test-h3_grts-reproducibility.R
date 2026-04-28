# test-h3_grts-reproducibility.R
# Tests that h3_grts() output is deterministic for a given seed and differs
# across seeds. The snapshot test in test-h3_grts-equivalence.R provides the
# stronger cross-platform stability anchor; these tests focus on within-session
# determinism.

test_that("same seed produces identical results across calls", {
  skip_if_not_installed("sf")
  wkb <- create_test_wkb()
  r1 <- h3_grts(wkb, n = 5, resolution = 4, seed = 42)
  r2 <- h3_grts(wkb, n = 5, resolution = 4, seed = 42)
  expect_equal(r1$cell,     r2$cell)
  expect_equal(r1$sort_key, r2$sort_key)
  expect_equal(r1$area_km2, r2$area_km2)
  expect_equal(r1$ip,       r2$ip)
})

test_that("different seeds produce different cell selections", {
  skip_if_not_installed("sf")
  wkb <- create_test_wkb()
  r1 <- h3_grts(wkb, n = 5, resolution = 4, seed = 1)
  r2 <- h3_grts(wkb, n = 5, resolution = 4, seed = 2)
  expect_false(identical(r1$cell, r2$cell))
})

test_that("area_correction = TRUE is reproducible with the same seed", {
  skip_if_not_installed("sf")
  wkb <- create_test_wkb()
  r1 <- h3_grts(wkb, n = 100, resolution = 6, area_correction = TRUE, seed = 7)
  r2 <- h3_grts(wkb, n = 100, resolution = 6, area_correction = TRUE, seed = 7)
  expect_equal(r1$cell, r2$cell)
  expect_equal(r1$sort_key, r2$sort_key)
  expect_equal(r1$ip, r2$ip)
})

test_that("area_correction = TRUE: different seeds produce different samples", {
  skip_if_not_installed("sf")
  wkb <- create_test_wkb()
  r1 <- h3_grts(wkb, n = 5, resolution = 4, area_correction = TRUE, seed = 7)
  r2 <- h3_grts(wkb, n = 5, resolution = 4, area_correction = TRUE, seed = 8)
  expect_false(identical(r1$cell, r2$cell))
})
