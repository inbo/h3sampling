#' Generate a GRTS sample on a continuous spatial domain
#'
#' @description
#' Generates a continuous Generalized Random Tessellation Stratified (GRTS)
#' sample from a `POLYGON` or `MULTIPOLYGON` without relying on discrete grid
#' approximations or rejection sampling.
#'
#' The algorithm generates a 1D spatial sequence (via interleaved 2D
#' Owen-scrambled Sobol sequences) that routes directly through the geometry's
#' exact area.
#' It traverses a lazy top-down quadtree in equal-area `(lon, sin(lat))` space,
#' dynamically clipping the polygon to compute exact intersection areas at every
#' branch.
#' This guarantees equal(-area) inclusion probabilities and (practically)
#' infinite continuous resolution.
#'
#' @param wkb Raw byte vector representing a Well-Known Binary (WKB) geometry
#'   in WGS84 (EPSG:4326). Must be a `POLYGON` or `MULTIPOLYGON`.
#' @param n Integer. The target number of sample locations to draw.
#' @param seed Numeric. A random seed passed to the Rust engine to ensure
#'   reproducible sampling. Must be a whole number `>= 0` and `< 2^53`.
#'   The seed controls the independent Owen scrambles applied to the spatial
#'   dimensions before Morton interleaving.
#'
#' @return A `data.frame` with `n` rows and three columns:
#' * `lon`         – Longitude of the sample point (WGS84 decimal degrees).
#' * `lat`         – Latitude  of the sample point (WGS84 decimal degrees).
#' * `sobol_index` – Index into the 1D GRTS stream at which this point was
#'                   generated. Useful for verifying prefix stability.
#'
#' The returned `data.frame` also carries the following attributes:
#' * `seed`            – The seed used (numeric).
#' * `n`               – Requested sample size (integer).
#' * `bbox`            – Named numeric vector `c(xmin, ymin, xmax, ymax)` of the
#'                       geometry bounding box in decimal degrees.
#' * `sobol_scanned`   – Because this is a zero-rejection algorithm, this will
#'                       always exactly equal `n`.
#' * `acceptance_rate` – Always `1.0`. The algorithm strictly routes points
#'                       into valid geometry, completely bypassing rejection.
#'
#' @section Exact Equal-Area Routing:
#' The algorithm maps the sequence to a quadtree rooted at global
#' longitude `[-180, 180]` and `sin(latitude)` `[-1, 1]`.
#' By operating in `sin(lat)` space, all cells are exactly equal-area on the
#' planetary sphere, eliminating polar distortion.
#' Polygon edges are automatically densified prior to transformation to limit
#' any spherical curvature containment error to < 1.4 m.
#'
#' @section Coordinate system & Antimeridian:
#' Input and output are in WGS84 (EPSG:4326) decimal degrees.
#' No external projection is required.
#' However, polygons straddling the antimeridian (180°/-180°) must be
#' pre-split into a `MULTIPOLYGON` (e.g., using `sf::st_wrap_dateline()`) before
#' conversion to WKB.
#'
#' @section Local Prefix Stability & Master Samples:
#' Calling with the same `seed` and a larger sample size `n2 > n1` on the
#' **same polygon** produces a sample whose first `n1` rows are identical to the
#' `n1` row result.
#'
#' Note: Because the 1D routing relies on the exact local area of the provided
#' polygon to achieve zero rejection, spatial consistency is **not** preserved
#' across overlapping but differently shaped polygons.
#' For globally consistent master samples, users should generate points over a
#' shared bounding box and spatially filter them in R.
#' Alternatively, use `bas_sobol()` for this purpose.
#'
#' @references
#' Robertson, B. L., et al. (2013). Spatial-balanced sampling in continuous
#' space. \emph{Environmental and Ecological Statistics}, 20(1), 149-167.
#'
#' Burley, B. (2020). Practical Hash-based Owen Scrambling.
#' \emph{Journal of Computer Graphics Techniques (JCGT)}, 9(4), 1-20.
#' \url{https://jcgt.org/published/0009/04/01/}
#'
#' Joe, S., & Kuo, F. Y. (2008). Constructing Sobol sequences with better
#' two-dimensional projections.
#' \emph{SIAM Journal on Scientific Computing}, 30(5), 2635-2654.
#' \doi{10.1137/070709359}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Load the built-in North Carolina dataset from the sf package
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#'
#' # Use the first county (Ashe) as our study area
#' study_area <- nc[1, ] |> st_geometry() |> st_as_binary(EWKB = FALSE)
#'
#' # Draw a spatially balanced sample of 50 points
#' sample_df <- grts_sobol(
#'   wkb  = study_area[[1]],
#'   n    = 50,
#'   seed = 2021
#' )
#'
#' # Prefix stability: extend to 100 points — first 50 rows are identical
#' sample_df_100 <- grts_sobol(
#'   wkb  = study_area[[1]],
#'   n    = 100,
#'   seed = 2021
#' )
#' stopifnot(identical(sample_df[1:50, c("lon", "lat")],
#'                     sample_df_100[1:50, c("lon", "lat")]))
#'
#' # Convert to sf and visualise
#' sample_sf <- st_as_sf(sample_df, coords = c("lon", "lat"), crs = 4326)
#' plot(st_geometry(nc[1, ]), main = "GRTS-Sobol Sample (Ashe County)")
#' plot(st_geometry(sample_sf), add = TRUE, pch = 20, col = "darkgreen")
#' }
grts_sobol <- function(wkb, n, seed) {
  # 0. Assertions
  stopifnot(
    "Input 'wkb' must be a raw byte vector representing WKB geometry." =
      is.raw(wkb),
    "`n` must be a whole number, > 0, not `NA`." =
      is_valid_n(n),
    "`seed` must be a whole number, >= 0, and < 2^53." =
      is_valid_seed(seed)
  )

  # Call the Rust engine
  res_df <- grts_sobol_sample_from_wkb(
    wkb_bytes = wkb,
    n = as.integer(n),
    seed = as.numeric(seed)
  )

  # Attach names to the bbox attribute
  bbox <- attr(res_df, "bbox")
  names(bbox) <- c("xmin", "ymin", "xmax", "ymax")
  attr(res_df, "bbox") <- bbox

  return(res_df)
}
