test_that("h3_grts output is strictly portable across OS architectures", {
  skip_if_not_installed("sf")
  skip_if_not_installed("h3o")
  nc_path <- system.file("shape/nc.shp", package = "sf")
  skip_if(nc_path == "", "sf built-in dataset not found")

  nc <- sf::st_read(nc_path, quiet = TRUE)
  ashe_wkb <- nc[1, ] |> sf::st_geometry() |> sf::st_as_binary()
  ashe_wkb <- ashe_wkb[[1]]

  # 1. Run a small, fixed sample
  target_n <- 5
  res_level <- 7
  test_seed <- 42

  sample_out <- h3_grts(
    wkb = ashe_wkb,
    n = target_n,
    resolution = res_level,
    containment = "centroid",
    seed = test_seed
  )

  # 2. The "Known Output"
  expected_strings <- c(
    "872a99369ffffff",
    "872a99254ffffff",
    "872a9922bffffff",
    "872a99361ffffff",
    "872a99249ffffff"
  )

  # 3. Strict assertion
  expect_equal(
    sample_out,
    expected_strings
  )
})
