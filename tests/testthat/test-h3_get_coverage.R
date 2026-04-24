# ---------------------------------------------------------------------------
# 1. Input Validation Tests
# ---------------------------------------------------------------------------
test_that("h3_get_coverage validates inputs correctly", {
  testthat::skip_if_not_installed("sf")
  valid_wkb <- create_test_wkb()

  # Invalid wkb
  expect_error(
    h3_get_coverage(wkb = "not a raw vector"),
    "Input 'wkb' must be a raw byte vector"
  )
  expect_error(
    h3_get_coverage(wkb = NULL),
    "Input 'wkb' must be a raw byte vector"
  )

  # Invalid resolution
  expect_error(h3_get_coverage(valid_wkb, resolution = -1))
  expect_error(h3_get_coverage(valid_wkb, resolution = 16))
  expect_error(h3_get_coverage(valid_wkb, resolution = "5"))

  # Invalid containment mode
  expect_error(
    h3_get_coverage(valid_wkb, containment = "inside"),
    "should be one of"
  )
})

# ---------------------------------------------------------------------------
# 2. Output Type & Bincode Format Tests
# ---------------------------------------------------------------------------
test_that("h3_get_coverage returns a raw bincode vector", {
  valid_wkb <- create_test_wkb()
  cells_raw <- h3_get_coverage(valid_wkb, resolution = 4)

  # The new bincode implementation returns raw bytes, not character strings
  expect_type(cells_raw, "raw")

  # Bincode serialization of a Vec<T> includes an 8-byte length prefix.
  # So even an empty coverage will be exactly 8 bytes.
  # Populated coverages will be > 8 bytes.
  expect_true(length(cells_raw) >= 8)
})

# ---------------------------------------------------------------------------
# 3. Logical/Algorithm Tests
# ---------------------------------------------------------------------------
test_that("h3_get_coverage respects resolution scaling", {
  valid_wkb <- create_test_wkb()

  # H3 grids scale logarithmically: higher resolution = more, smaller cells
  res3_bytes <- h3_get_coverage(valid_wkb, resolution = 3)
  res4_bytes <- h3_get_coverage(valid_wkb, resolution = 4)

  # A higher resolution coverage should serialize to a larger raw byte vector
  expect_true(length(res4_bytes) > length(res3_bytes))
})

test_that("h3_get_coverage accepts all valid containment modes", {
  valid_wkb <- create_test_wkb()

  # Ensure the Rust parsing logic successfully accepts all 4 modes
  expect_no_error(h3_get_coverage(valid_wkb, containment = "centroid"))
  expect_no_error(h3_get_coverage(valid_wkb, containment = "intersect"))
  expect_no_error(h3_get_coverage(valid_wkb, containment = "boundary"))
  expect_no_error(h3_get_coverage(valid_wkb, containment = "covers"))
})

# ---------------------------------------------------------------------------
# 4. Snapshot / Stability Tests
# ---------------------------------------------------------------------------
test_that("h3_get_coverage produces stable byte payloads across builds", {
  fixed_wkb <- create_test_wkb()

  # Generate the payload
  cells_raw <- h3_get_coverage(
    fixed_wkb, resolution = 4, containment = "centroid"
  )

  # style = "serialize" is perfect for raw bytes
  expect_snapshot_value(cells_raw, style = "serialize")
})
