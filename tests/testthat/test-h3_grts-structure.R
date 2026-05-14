# test-h3_grts-structure.R
# Tests for the structure, column types, and attributes of h3_grts() output.
# These do not depend on specific cell values — they hold for any valid input.

# ---------------------------------------------------------------------------
# 1. Return type and columns
# ---------------------------------------------------------------------------

test_that("h3_grts returns a data.frame with correct columns", {
  skip_if_not_installed("sf")
  result <- h3_grts(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_s3_class(result, "data.frame")
  expect_equal(
    names(result),
    c("cell", "sort_key", "effective_key", "area_m2", "ip")
  )
})

test_that("h3_grts columns have correct types", {
  skip_if_not_installed("sf")
  result <- h3_grts(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_type(result$cell, "character")
  expect_type(result$sort_key, "character")
  expect_type(result$area_m2, "double")
  expect_type(result$ip, "double")
})

test_that("cell IDs are valid lowercase hex without zero-padding", {
  skip_if_not_installed("sf")
  result <- h3_grts(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_true(all(grepl("^[0-9a-f]+$", result$cell)))
  # Canonical H3 format: no leading zeros (unless the value itself is zero)
  expect_false(any(grepl("^0[0-9a-f]", result$cell)))
})

test_that("sort_key is zero-padded to exactly 16 lowercase hex characters", {
  skip_if_not_installed("sf")
  result <- h3_grts(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_true(all(nchar(result$sort_key) == 16L))
  expect_true(all(grepl("^[0-9a-f]{16}$", result$sort_key)))
})

test_that("rows are ordered by ascending sort_key", {
  skip_if_not_installed("sf")
  # Zero-padded hex: lexicographic order == numeric order of the u64 sort key
  result <- h3_grts(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_equal(result$sort_key, sort(result$sort_key))
})

test_that("cell IDs are unique within a sample", {
  skip_if_not_installed("sf")
  result <- h3_grts(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_equal(length(unique(result$cell)), nrow(result))
})

test_that("sort keys are unique within a sample", {
  skip_if_not_installed("sf")
  result <- h3_grts(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_equal(length(unique(result$sort_key)), nrow(result))
})

# ---------------------------------------------------------------------------
# 2. Attributes
# ---------------------------------------------------------------------------

test_that("all design attributes are present with correct types", {
  skip_if_not_installed("sf")
  result <- h3_grts(create_test_wkb(), n = 5, resolution = 4,
                    containment = "centroid", seed = 42)
  expect_type(attr(result, "n_cells"), "integer")
  expect_type(attr(result, "sum_area_m2"), "double")
  expect_type(attr(result, "n"), "integer")
  expect_type(attr(result, "resolution"), "integer")
  expect_type(attr(result, "containment"), "character")
  expect_type(attr(result, "area_correction"), "logical")
})

test_that("attribute n equals nrow(result)", {
  skip_if_not_installed("sf")
  result <- h3_grts(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_equal(attr(result, "n"), nrow(result))
})

test_that("attribute containment echoes the input argument", {
  skip_if_not_installed("sf")
  result <- h3_grts(create_test_wkb(), n = 5, resolution = 4,
                    containment = "centroid", seed = 42)
  expect_equal(attr(result, "containment"), "centroid")
})

test_that("area_correction attribute is FALSE by default", {
  skip_if_not_installed("sf")
  result <- h3_grts(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_false(attr(result, "area_correction"))
})

test_that("n_cells >= n", {
  skip_if_not_installed("sf")
  result <- h3_grts(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_gte(attr(result, "n_cells"), 5L)
})

test_that("sum_area_m2 is positive", {
  skip_if_not_installed("sf")
  result <- h3_grts(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_gt(attr(result, "sum_area_m2"), 0)
})
