# test-reproducibility.R
# Tests that sampler output is deterministic for a given seed, differs
# across seeds and seeds are cross-platform stable

# ---------------------------------------------------------------------------
# 1. Snapshot — cross-OS portability anchor
# ---------------------------------------------------------------------------

test_that("output is strictly portable across OS architectures", {
  skip_if_not_installed("sf")
  test_wkb <- create_test_wkb()
  h3_grts_out <- h3_grts(
    wkb         = test_wkb,
    n           = 5,
    resolution  = 4,
    containment = "centroid",
    seed        = 42
  )
  # style = "serialize" captures exact binary equivalence including column
  # types, float bit patterns, and attribute values.
  expect_snapshot_value(h3_grts_out, style = "serialize")

  grts_sobol_out <- grts_sobol(
    wkb         = test_wkb,
    n           = 5,
    seed        = 42
  )
  expect_snapshot_value(grts_sobol_out, style = "serialize")

  bas_sobol_out <- bas_sobol(
    wkb         = test_wkb,
    n           = 5,
    seed        = 42
  )
  expect_snapshot_value(bas_sobol_out, style = "serialize")
})

# ---------------------------------------------------------------------------
# 2. Within session test
# ---------------------------------------------------------------------------

test_that("same seed produces identical results across calls", {
  skip_if_not_installed("sf")
  wkb <- create_test_wkb()
  r1 <- h3_grts(wkb, n = 5, resolution = 4, seed = 42)
  r2 <- h3_grts(wkb, n = 5, resolution = 4, seed = 42)
  expect_equal(r1, r2)

  r1 <- h3_grts(wkb, n = 100, resolution = 6, area_correction = TRUE, seed = 7)
  r2 <- h3_grts(wkb, n = 100, resolution = 6, area_correction = TRUE, seed = 7)
  expect_equal(r1, r2)

  r1 <- grts_sobol(wkb, n = 5, seed = 42)
  r2 <- grts_sobol(wkb, n = 5, seed = 42)
  expect_equal(r1, r2)

  r1 <- bas_sobol(wkb, n = 5, seed = 42)
  r2 <- bas_sobol(wkb, n = 5, seed = 42)
  expect_equal(r1, r2)
})

test_that("different seeds produce different selections", {
  skip_if_not_installed("sf")
  wkb <- create_test_wkb()
  r1 <- h3_grts(wkb, n = 5, resolution = 4, seed = 1)
  r2 <- h3_grts(wkb, n = 5, resolution = 4, seed = 2)
  expect_false(identical(r1$cell, r2$cell))

  r1 <- h3_grts(wkb, n = 5, resolution = 4, area_correction = TRUE, seed = 7)
  r2 <- h3_grts(wkb, n = 5, resolution = 4, area_correction = TRUE, seed = 8)
  expect_false(identical(r1$cell, r2$cell))

  r1 <- grts_sobol(wkb, n = 5, seed = 1)
  r2 <- grts_sobol(wkb, n = 5, seed = 2)
  expect_false(identical(r1$lon, r2$lon))

  r1 <- bas_sobol(wkb, n = 5, seed = 1)
  r2 <- bas_sobol(wkb, n = 5, seed = 2)
  expect_false(identical(r1$lon, r2$lon))
})
