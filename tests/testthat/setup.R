# ---------------------------------------------------------------------------
# Setup helper: Create a valid WKB polygon for testing
# ---------------------------------------------------------------------------
create_test_wkb <- function() {
  # A simple 1x1 degree square near the equator
  poly <- sf::st_polygon(
    list(matrix(c(0,0, 0,1, 1,1, 1,0, 0,0), ncol=2, byrow=TRUE))
  )
  poly_sfc <- sf::st_sfc(poly, crs = 4326)
  sf::st_as_binary(poly_sfc, EWKB = FALSE)[[1]]
}
