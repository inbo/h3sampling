# test-h3_grts-equivalence.R
# Tests for cross-platform stability (snapshot) and equivalence between the
# WKB and pre-computed cells input pathways.

# ---------------------------------------------------------------------------
# 1. Snapshot — cross-OS portability anchor
# ---------------------------------------------------------------------------

test_that("h3_grts output is strictly portable across OS architectures", {
  skip_if_not_installed("sf")
  sample_out <- h3_grts(
    wkb         = create_test_wkb(),
    n           = 5,
    resolution  = 4,
    containment = "centroid",
    seed        = 42
  )
  # style = "serialize" captures exact binary equivalence including column
  # types, float bit patterns, and attribute values.
  expect_snapshot_value(sample_out, style = "serialize")
})

# ---------------------------------------------------------------------------
# 2. WKB vs pre-computed cells pathway equivalence
# ---------------------------------------------------------------------------

test_that("WKB and cells pathways produce identical samples", {
  skip_if_not_installed("sf")
  wkb <- create_test_wkb()
  n <- 10
  res <- 5
  seed <- 123

  out_wkb <- h3_grts(
    wkb         = wkb,
    n           = n,
    resolution  = res,
    containment = "centroid",
    seed        = seed
  )

  coverage_bytes <- h3_get_coverage(
    wkb         = wkb,
    resolution  = res,
    containment = "centroid"
  )
  out_cells <- h3_grts(
    cells = coverage_bytes,
    n     = n,
    seed  = seed
  )

  # The cells path records containment as "pre-computed" — nullify before
  # comparing so the structural equality check is not confused by this
  # expected metadata difference.
  expect_equal(attr(out_wkb,   "containment"), "centroid")
  expect_equal(attr(out_cells, "containment"), "pre-computed")
  attr(out_wkb,   "containment") <- NULL
  attr(out_cells, "containment") <- NULL

  expect_equal(out_wkb, out_cells)
})

test_that("cells pathway: sort keys for shared cells match the WKB pathway", {
  skip_if_not_installed("sf")
  # The sort key for any H3 cell is a pure function of (cell, seed).
  # A cell appearing in both samples must carry the same sort_key.
  wkb  <- create_test_wkb()
  seed <- 42

  out_wkb <- h3_grts(
    wkb, n = 5, resolution = 4,
    containment = "centroid", seed = seed
  )
  coverage_bytes <- h3_get_coverage(
    wkb, resolution = 4, containment = "centroid"
  )
  out_cells <- h3_grts(cells = coverage_bytes, n = 5, seed = seed)

  shared <- intersect(out_wkb$cell, out_cells$cell)

  keys_wkb   <- out_wkb$sort_key[match(shared, out_wkb$cell)]
  keys_cells <- out_cells$sort_key[match(shared, out_cells$cell)]
  expect_equal(keys_wkb, keys_cells)
})

# ---------------------------------------------------------------------------
# 3. Input validation (dispatch and type checks on R side)
# ---------------------------------------------------------------------------

test_that("h3_grts errors when neither wkb nor cells is provided", {
  expect_error(
    h3_grts(n = 5),
    "Either 'wkb' or 'cells' must be provided"
  )
})

test_that("h3_grts errors when cells is a character vector instead of raw", {
  expect_error(
    h3_grts(cells = c("842a993ffffffff"), n = 5),
    "must be a raw byte vector"
  )
})

test_that("h3_grts errors for invalid resolution", {
  skip_if_not_installed("sf")
  expect_error(h3_grts(create_test_wkb(), n = 3, resolution = -1))
  expect_error(h3_grts(create_test_wkb(), n = 3, resolution = 16))
})

test_that("h3_grts errors for invalid containment mode", {
  skip_if_not_installed("sf")
  expect_error(
    h3_grts(create_test_wkb(), n = 3, containment = "overlap"),
    "should be one of"
  )
})

test_that("h3_grts containment ordering: intersect >= centroid cells", {
  skip_if_not_installed("sf")
  wkb <- create_test_wkb()
  n_intersect <- attr(
    h3_grts(wkb, n = 1, resolution = 4, containment = "intersect"),
    "n_cells"
  )
  n_centroid <- attr(
    h3_grts(wkb, n = 1, resolution = 4, containment = "centroid"),
    "n_cells"
  )
  expect_gte(n_intersect, n_centroid)
})
