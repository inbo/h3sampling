#' Generate a Spatially Balanced GRTS Sample using H3 (Equal Probabilities)
#'
#' @description
#' Generates a continuous Generalized Random-Tessellation Stratified (GRTS)
#' sample.
#' This function tessellates a study area polygon into an H3 hexagonal grid and
#' orders the cells using a Reverse-Radix master sample key.
#' It assumes equal inclusion probabilities and draws the first `n` elements.
#'
#' @param resolution The H3 resolution to use (0-15).
#' @param containment The containment mode to use. One of "centroid",
#' "intersect", "boundary", or "covers".
#' @param wkb Raw byte vector representing Well-Known Binary (WKB) geometry in
#' WGS84 (EPSG:4326).
#' @param area_correction If `FALSE` (default), all H3 cells have equal
#'  inclusion probability pi_i = n / N.
#'  Cell areas vary by up to ~1.9x across the icosahedral projection;
#'  this is ignored in the sampling step but can be corrected at the estimation
#'  stage using the returned `area_km2` and `ip` columns.
#'  If `TRUE`, rejection sampling is applied so that each cell is accepted with
#'  probability area_i / max_area, making the effective inclusion probability
#'  proportional to area: pi_i = n * area_i / sum(area).
#'  A pre-pass over the coverage is required to find the maximum cell area, so
#'  this mode incurs roughly 2x the runtime of the default mode.
#' @param n Integer. The target number of sample locations to draw.
#' @param seed Numeric. A random seed passed to the Rust engine to ensure a
#'   reproducible hierarchical shuffle.
#'
#'
#' @return A `data.frame` with four columns and `n` rows, ordered by ascending
#'  sort key (GRTS visiting order):
#' * `cell`     – H3 cell identifier as a lowercase hexadecimal string.
#' * `sort_key` – GRTS sort key as a zero-padded 16-character lowercase hex
#'                string. Zero-padding ensures correct lexicographic ordering
#'                (e.g. for merging tile results).
#' * `area_m2` – True area of the cell in m².
#' * `ip`      – Inclusion probability pi_i of the sampled cell:
#'                  area_correction = FALSE: pi_i = n / N  (constant)
#'                  area_correction = TRUE:  pi_i = n * area_i / sum(area)
#' Use `h3o::h3_to_points(h3o::h3_from_strings())` to convert the `cell` strings
#' back to `sf` geometries.
#' The returned `data.frame` also carries the following attributes which
#' provide the quantities needed to reconstruct design-based estimators:
#' * `n_cells`     – Total number of H3 cells N covering the study area.
#'                   Used to compute pi_i = n / N.
#' * `sum_area_m2`– Sum of all cell areas in m² (i.e. the approximate
#'                   area of the study area as seen by the H3 grid).
#' * `n`           – Requested sample size.
#' * `resolution`  – H3 resolution used.
#' * `containment` – Containment mode used.
#' * `area_correction` – Whether area-proportional sampling was applied.
#' @export
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(h3o)
#'
#' # Load the built-in North Carolina dataset from the sf package
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#'
#' # Use the first county (Ashe) as our study area
#' study_area <- nc[1, ] |> st_geometry() |> st_as_binary(EWKB = FALSE)
#'
#' # Draw a sample of 25 points at resolution 7
#' sample_df <- h3_grts(
#'   wkb = study_area[[1]],
#'   n = 25,
#'   resolution = 7,
#'   containment = "centroid",
#'   seed = 123
#' )
#'
#' # Convert the H3 indices back into spatial points
#' sample_df$geom <- h3_to_points(h3_from_strings(sample_df$cell))
#' sample_df <- st_as_sf(sample_df)
#'
#' # Visualize the spatially balanced sample
#' plot(st_geometry(nc[1,]), main = "H3-GRTS Sample (Ashe County)")
#' plot(st_geometry(sample_df), add = TRUE, pch = 20, col = "red")
#' }
h3_grts <- function(
  wkb,
  n,
  resolution = 5,
  containment = c("centroid", "intersect", "boundary", "covers"),
  seed = 42,
  area_correction = FALSE
) {
  # 0. Assertions
  if (!is.raw(wkb)) {
    stop(
      "Input 'wkb' must be a raw byte vector representing Well-Known Binary (WKB) geometry.\n", # nolint
      "Hint: If you have an sf polygon 'poly', you can generate this using:\n",
      "  wkb <- sf::st_as_binary(sf::st_combine(poly), EWKB = FALSE)[[1]]"
    )
  }
  stopifnot(
    is.numeric(n),
    is.numeric(resolution),
    is.numeric(seed),
    is.logical(area_correction) & !is.na(area_correction),
    n > 0,
    seed > 0 & seed < 2^53,
    resolution >= 0 & resolution <= 15
  )
  containment <- match.arg(containment)

  # 1. Call the Rust Pipeline
  # Rust handles H3 tessellation and GRTS sorting internally
  selected_cells_str <- generate_grts_sample(
    wkb_bytes = wkb,
    res = as.integer(resolution),
    containment = containment,
    n = as.integer(n),
    global_seed = as.numeric(seed),
    area_correction = area_correction
  )

  return(selected_cells_str)
}
