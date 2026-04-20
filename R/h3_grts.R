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
#' @param n Integer. The target number of sample locations to draw.
#' @param seed Numeric. A random seed passed to the Rust engine to ensure a
#'   reproducible hierarchical shuffle.
#'
#'
#' @return A data.frame with `n` rows and two columns.
#'   The first column contains selected H3 cell indices as a character vector.
#'   The second column contains the sort index.
#'   Use `h3o::h3_to_points(h3o::h3_from_strings())` to convert H3 cell indices
#'   back to `sf` geometries.
#'
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
#' sample_cells <- h3_grts(
#'   wkb = study_area[[1]],
#'   n = 25,
#'   resolution = 7,
#'   containment = "centroid",
#'   seed = 123
#' )
#'
#' # Convert the H3 indices back into spatial points
#' sample_points <- h3_to_points(h3_from_strings(sample_cells))
#'
#' # Visualize the spatially balanced sample
#' plot(st_geometry(nc[1,]), main = "H3-GRTS Sample (Ashe County)")
#' plot(st_geometry(sample_points), add = TRUE, pch = 20, col = "red")
#' }
h3_grts <- function(
    wkb,
    n,
    resolution = 5,
    containment = c("centroid", "intersect", "boundary", "covers"),
    seed = 42
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
    n > 0,
    seed > 0,
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
    global_seed = as.numeric(seed)
  )

  return(selected_cells_str)
}
