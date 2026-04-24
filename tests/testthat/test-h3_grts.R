# ---------------------------------------------------------------------------
# 1. Snapshot / Portability Test
# ---------------------------------------------------------------------------
test_that("h3_grts output is strictly portable across OS architectures", {
  test_wkb <- create_test_wkb()

  sample_out <- h3_grts(
    wkb = test_wkb,
    n = 5,
    resolution = 4,
    containment = "centroid",
    seed = 42
  )

  # Using style = "serialize" guarantees exact binary equivalence
  # across Mac/Linux/Windows,
  # including column types, attributes, and bit-level float precision.
  expect_snapshot_value(sample_out, style = "serialize")
})

# ---------------------------------------------------------------------------
# 2. Pathway Equivalence Test (WKB vs. Pre-computed Cells)
# ---------------------------------------------------------------------------
test_that("WKB and Cells pathways produce identical samples", {
  test_wkb <- create_test_wkb()
  target_n <- 10
  res <- 5
  test_seed <- 123

  # 1. Run the standard WKB on-the-fly pathway
  out_wkb <- h3_grts(
    wkb = test_wkb,
    n = target_n,
    resolution = res,
    containment = "centroid",
    seed = test_seed
  )

  # 2. Run the pre-computed cells pathway
  coverage_bytes <- h3_get_coverage(
    wkb = test_wkb,
    resolution = res,
    containment = "centroid"
  )
  out_cells <- h3_grts(
    cells = coverage_bytes,
    n = target_n,
    seed = test_seed
  )

  # 3. Handle expected attribute differences
  # The WKB path uses "centroid", but the cells path hardcodes "pre-computed"
  expect_equal(attr(out_wkb, "containment"), "centroid")
  expect_equal(attr(out_cells, "containment"), "pre-computed")

  # Nullify the containment attribute so we can strictly compare the rest
  attr(out_wkb, "containment") <- NULL
  attr(out_cells, "containment") <- NULL

  # 4. Strict assertion
  expect_equal(out_wkb, out_cells)
})

# ---------------------------------------------------------------------------
# 3. Input Validation & Edge Cases
# ---------------------------------------------------------------------------
test_that("h3_grts validates inputs and handles missing args correctly", {
  test_wkb <- create_test_wkb()

  # Missing both geometry and cells
  expect_error(
    h3_grts(n = 5),
    "Either 'wkb' or 'cells' must be provided"
  )

  # Invalid sample size
  expect_error(
    h3_grts(wkb = test_wkb, n = -5),
    "n > 0 is not TRUE"
  )

  # Invalid cell vector type (expects raw bytes, not character)
  expect_error(
    h3_grts(cells = c("842a993ffffffff"), n = 5),
    "must be a raw byte vector"
  )

  # Target n exceeds the number of available cells in the polygon
  expect_error(
    h3_grts(wkb = test_wkb, n = 50000, resolution = 2),
    "panicked"
  )
})

# ---------------------------------------------------------------------------
# 4. Area Correction Execution Test
# ---------------------------------------------------------------------------
test_that("h3_grts executes successfully with area_correction enabled", {
  test_wkb <- create_test_wkb()

  out <- h3_grts(
    wkb = test_wkb,
    n = 10,
    resolution = 5,
    area_correction = TRUE
  )

  expect_s3_class(out, "data.frame")
  expect_true(attr(out, "area_correction"))
  # Inclusion probabilities should not be universally equal
  expect_true(is.numeric(out$ip))
  expect_all_false(out$ip == mean(out$ip))
})
