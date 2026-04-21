use extendr_api::prelude::*;
use h3o::{CellIndex, Resolution};
use h3o::geom::{TilerBuilder, ContainmentMode};
use geo::Geometry;
use geozero::wkb::Wkb;
use geozero::ToGeo;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_pcg::Pcg64Mcg;
use std::collections::{BinaryHeap, HashMap};

// ---------------------------------------------------------------------------
// H3 bit-layout constants (from the H3 spec)
// ---------------------------------------------------------------------------

/// Bit position where the base cell index starts in an H3 index.
const H3_BASE_CELL_SHIFT: u32 = 45;

/// 7-bit mask to extract the base cell number (values 0–121).
const H3_BASE_CELL_MASK: u64 = 0b111_1111;

/// Maximum valid base cell index.
const H3_MAX_BASE_CELL: usize = 121;

/// Number of base cells in the H3 grid.
const H3_NUM_BASE_CELLS: usize = 122;

/// 3-bit mask to extract a single resolution digit (child index 0–6).
const H3_DIGIT_MASK: u64 = 0b111;

/// Number of bits used to store the base cell in the sort key.
const SORT_KEY_BASE_BITS: u32 = 7;

/// Number of bits used per resolution level in the sort key.
const SORT_KEY_BITS_PER_RES: u32 = 3;

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Generate a GRTS (Generalized Random Tessellation Stratified) sample of H3
/// cell indices covering a polygon.
///
/// # Arguments
/// * `wkb_bytes`  – Study area geometry encoded as Well-Known Binary.
/// * `res`        – H3 resolution (0 = coarsest, 15 = finest).
/// * `containment`– One of `"intersect"`, `"centroid"`, `"boundary"`, `"covers"`.
/// * `n`          – Desired sample size.
/// * `global_seed`– Reproducibility seed (passed as f64 from R; values above
///                  2^53 lose precision — prefer seeds within that range).
///
/// # Returns
/// A `data.frame` with three columns, ordered by ascending sort key (GRTS
/// visiting order):
/// * `cell`     – H3 cell identifier as a lowercase hexadecimal string.
/// * `sort_key` – GRTS sort key as a zero-padded 16-character lowercase hex
///                string. Zero-padding ensures correct lexicographic ordering
///                on the R side (e.g. for merging tile results).
/// * `area_m2` – True area of the cell in m². Because GRTS samples cells
///                with equal inclusion probability (π_i = n / N), area is NOT
///                part of the sampling mechanism. It enters only at the
///                estimation stage: to estimate a population total T, scale
///                each observation y_i by its cell area A_i and divide by π_i:
///                  T_hat = Σ (y_i * A_i) / π_i
///                For the Hájek mean estimator the π_i terms cancel and the
///                result is simply an area-weighted average of the y_i values.
///
/// The returned `data.frame` also carries the following attributes which
/// provide the quantities needed to reconstruct any design-based estimator:
/// * `n_cells`     – Total number of H3 cells N covering the study area.
///                   Used to compute π_i = n / N.
/// * `sum_area_m2`– Sum of all cell areas in m² (i.e. the approximate
///                   area of the study area as seen by the H3 grid).
/// * `n`           – Requested sample size.
/// * `resolution`  – H3 resolution used.
/// * `containment` – Containment mode used.
#[extendr]
fn generate_grts_sample(
    wkb_bytes: &[u8],
    res: i32,
    containment: &str,
    n: i32,
    global_seed: f64,
) -> extendr_api::Result<Robj> {
    // NOTE: f64 → u64 truncates; seeds above 2^53 lose precision.
    let seed_u64 = global_seed as u64;
    let target_n = n as usize;

    // 1. Decode WKB and set up the H3 tiler --------------------------------
    let geometry = decode_wkb(wkb_bytes)?;
    let containment_mode = parse_containment(containment)?;
    let resolution = Resolution::try_from(res as u8)
        .map_err(|_| Error::Other(format!("Invalid H3 resolution: {}. Must be 0-15.", res)))?;

    let mut tiler = TilerBuilder::new(resolution)
        .containment_mode(containment_mode)
        .build();

    feed_geometry(&mut tiler, geometry)?;

    // 2. Shuffle base cells (top level of the GRTS address) ----------------
    // This randomises which continental regions receive low sort keys,
    // ensuring global spatial balance across the entire Earth surface.
    let base_cells = shuffled_base_cells(seed_u64);

    // 3. Stream cells through the heap reservoir ----------------------------
    // We use a max-heap of capacity `n` so we never materialise the full
    // coverage. Memory usage is O(n) regardless of total cell count.
    // Time complexity: O(total_cells x log n).
    //
    // The heap stores (sort_key, h3_index, area_bits) tuples. BinaryHeap is a
    // max-heap, so peek() always returns the *largest* key — exactly what we
    // want to evict when a smaller key arrives. area_bits is the cell area
    // (m²) bit-cast to u64 so it can be stored without a separate Vec.
    let mut heap: BinaryHeap<(u64, u64, u64)> = BinaryHeap::with_capacity(target_n + 1);

    // Cache permutations keyed by parent_id so siblings share a single RNG
    // initialisation instead of each independently re-seeding.
    let mut perm_cache: HashMap<u64, [u64; 7]> = HashMap::new();

    let mut total_cells: usize = 0;
    // Accumulate the total area across ALL cells in the coverage. This is
    // stored as the `sum_area_m2` attribute and is the denominator for
    // area-weighted estimators: μ_Hajek = Σ(y_i * A_i) / Σ A_i.
    let mut sum_all_areas: f64 = 0.0;

    for cell in tiler.into_coverage() {
        total_cells += 1;
        let cell_area = cell.area_m2();
        sum_all_areas += cell_area;
        let sort_key = compute_sort_key(cell, res, seed_u64, &base_cells, &mut perm_cache)?;
        reservoir_push(&mut heap, target_n, sort_key, u64::from(cell), cell_area);
    }

    // 4. Validate sample size -----------------------------------------------
    if total_cells < target_n {
        return Err(Error::Other(format!(
            "Target n ({}) exceeds the number of available cells ({}). \
             Try a higher resolution or a larger study area.",
            target_n, total_cells
        )));
    }

    // 5. Extract results in GRTS order --------------------------------------
    // BinaryHeap::into_sorted_vec() drains the max-heap in ascending order,
    // which corresponds directly to the GRTS visiting sequence.
    //
    // Sort keys are formatted as zero-padded 16-character hex strings so that
    // lexicographic ordering on the R side is identical to numeric ordering of
    // the underlying u64 values. This is essential for correct cross-tile
    // merging: order(sort_key) in R will give the true GRTS sequence.
    let sorted: Vec<(u64, u64, u64)> = heap.into_sorted_vec();

    let cells: Vec<String> = sorted
        .iter()
        .map(|&(_, idx, _)| format!("{:x}", idx))
        .collect();

    let sort_keys: Vec<String> = sorted
        .iter()
        .map(|&(key, _, _)| format!("{:016x}", key))
        .collect();

    // Extract per-cell areas for the sampled cells. Area enters only at the
    // estimation stage — it is NOT part of the equal-probability sampling
    // mechanism. π_i = n / N is constant for all cells.
    let areas_m2: Vec<f64> = sorted
        .iter()
        .map(|&(_, _, area_bits)| f64::from_bits(area_bits))
        .collect();

    // 6. Build the data.frame ----------------------------------------------
    let df = data_frame!(
        cell     = cells,
        sort_key = sort_keys,
        area_m2 = areas_m2
    );

    // 7. Attach design attributes -------------------------------------------
    // These provide everything needed to construct any design-based estimator
    // without baking assumptions into the returned object:
    //
    //   π_i  = n / n_cells           (equal for all sampled cells)
    //   T_hat = Σ (y_i * area_m2_i) / π_i
    //   μ_Hajek = Σ (y_i * area_m2_i / π_i) / Σ (area_m2_i / π_i)
    //           = Σ (y_i * area_m2_i) / Σ area_m2_i   (π_i cancels)
    let mut df_robj: Robj = df.into();
    df_robj.set_attrib("n_cells",      total_cells as i32)?;
    df_robj.set_attrib("sum_area_m2", sum_all_areas)?;
    df_robj.set_attrib("n",            n)?;
    df_robj.set_attrib("resolution",   res)?;
    df_robj.set_attrib("containment",  containment)?;

    Ok(df_robj)
}

// ---------------------------------------------------------------------------
// Sort key computation
// ---------------------------------------------------------------------------

/// Compute the GRTS sort key for a single H3 cell.
///
/// The sort key is a packed u64 built top-down through the H3 hierarchy:
///
/// ```text
/// bits 63..52  unused
/// bits 51..45  shuffled base cell  (7 bits, SORT_KEY_BASE_BITS)
/// bits 44..42  shuffled digit at res 1  (3 bits each)
/// bits 41..39  shuffled digit at res 2
/// ...
/// bits (44 - 3*(res-1))..(42 - 3*(res-1))  shuffled digit at target res
/// ```
///
/// Each digit is the child index (0–6) within its parent, remapped through a
/// per-parent pseudo-random permutation. Pentagon parents only have 6 valid
/// children (digits 0–5), so a 6-element permutation is used for those.
fn compute_sort_key(
    cell: CellIndex,
    res: i32,
    seed_u64: u64,
    base_cells: &[u64; H3_NUM_BASE_CELLS],
    perm_cache: &mut HashMap<u64, [u64; 7]>,
) -> extendr_api::Result<u64> {
    let h3_index = u64::from(cell);

    // -- Base cell (most significant part of the key) ----------------------
    let base_cell_raw = ((h3_index >> H3_BASE_CELL_SHIFT) & H3_BASE_CELL_MASK) as usize;
    let remapped_base = if base_cell_raw <= H3_MAX_BASE_CELL {
        base_cells[base_cell_raw]
    } else {
        base_cell_raw as u64
    };
    let mut sort_key: u64 = remapped_base;

    // -- Resolution digits (one per level from 1 to res) -------------------
    for r in 1..=res {
        // Use the h3o API to get the parent cell at resolution r-1, which is
        // more robust than manual bit-masking.
        let parent_res = Resolution::try_from((r - 1) as u8)
            .map_err(|_| Error::Other(format!("Invalid parent resolution: {}", r - 1)))?;
        let parent_cell = cell
            .parent(parent_res)
            .ok_or_else(|| Error::Other(format!("Could not get parent at res {}", r - 1)))?;
        let parent_id = u64::from(parent_cell);

        // Retrieve (or compute and cache) the permutation for this parent.
        // Pentagon parents have 6 children; all others have 7.
        let perm = perm_cache.entry(parent_id).or_insert_with(|| {
            let mut rng = Pcg64Mcg::seed_from_u64(parent_id.wrapping_add(seed_u64));
            let mut p = [0u64, 1, 2, 3, 4, 5, 6];
            // For pentagons the digit 1 (the "deleted" direction) is unused,
            // but we still shuffle all 7 slots — the deleted digit will simply
            // never be observed in practice.
            p.shuffle(&mut rng);
            p
        });

        // Extract the raw child digit for resolution r from the H3 index.
        // Each resolution occupies 3 bits; the most significant digit sits at
        // bit offset 44 (= H3_BASE_CELL_SHIFT - 1 - 2).
        let h3_offset = H3_BASE_CELL_SHIFT
            .checked_sub(3 * r as u32)
            .ok_or_else(|| Error::Other("Resolution digit offset underflow".into()))?;
        let digit = ((h3_index >> h3_offset) & H3_DIGIT_MASK) as usize;

        if digit < 7 {
            let remapped_digit = perm[digit];
            // Pack the remapped digit into the sort key.
            // Layout: 7 base-cell bits, then 3 bits per resolution level.
            let key_offset = SORT_KEY_BASE_BITS + SORT_KEY_BITS_PER_RES * (r as u32 - 1);
            sort_key |= remapped_digit << key_offset;
        }
    }

    // Sanity check: at res 15 we need 7 + 15*3 = 52 bits, well within u64.
    debug_assert!(
        SORT_KEY_BASE_BITS + SORT_KEY_BITS_PER_RES * res as u32 <= 64,
        "Sort key overflows u64 at resolution {}",
        res
    );

    Ok(sort_key)
}

// ---------------------------------------------------------------------------
// Heap reservoir helper
// ---------------------------------------------------------------------------

/// Push `(sort_key, h3_index, area_bits)` into the reservoir heap, evicting
/// the maximum element if the heap is already full and the new key is smaller.
/// Cell area is carried through so weights can be computed after the loop
/// without re-querying h3o.
#[inline]
fn reservoir_push(
    heap: &mut BinaryHeap<(u64, u64, u64)>,
    target_n: usize,
    sort_key: u64,
    h3_index: u64,
    cell_area: f64,
) {
    let area_bits = cell_area.to_bits();
    if heap.len() < target_n {
        heap.push((sort_key, h3_index, area_bits));
    } else if let Some(&(max_key, _, _)) = heap.peek() {
        if sort_key < max_key {
            heap.push((sort_key, h3_index, area_bits));
            heap.pop();
        }
    }
}

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

/// Decode a WKB byte slice into a `geo::Geometry`.
fn decode_wkb(wkb_bytes: &[u8]) -> extendr_api::Result<Geometry> {
    Wkb(wkb_bytes.to_vec())
        .to_geo()
        .map_err(|e| Error::Other(format!("WKB parse failed: {}", e)))
}

/// Parse the containment mode string into an `h3o` `ContainmentMode`.
fn parse_containment(containment: &str) -> extendr_api::Result<ContainmentMode> {
    match containment {
        "intersect" => Ok(ContainmentMode::IntersectsBoundary),
        "centroid"  => Ok(ContainmentMode::ContainsCentroid),
        "boundary"  => Ok(ContainmentMode::ContainsBoundary),
        "covers"    => Ok(ContainmentMode::Covers),
        other => Err(Error::Other(format!(
            "Unknown containment mode: '{}'. \
             Expected one of: intersect, centroid, boundary, covers.",
            other
        ))),
    }
}

/// Feed a `geo::Geometry` into an h3o tiler, recursing into collections.
/// Returns an error (rather than silently skipping) if unexpected geometry
/// types such as LineStrings are encountered, as these likely indicate an
/// upstream data pipeline issue.
fn feed_geometry(
    tiler: &mut h3o::geom::Tiler,
    geometry: Geometry,
) -> extendr_api::Result<()> {
    match geometry {
        Geometry::Polygon(p) => {
            tiler.add(p).map_err(|e| Error::Other(format!("Invalid polygon: {}", e)))?;
        }
        Geometry::MultiPolygon(mp) => {
            for p in mp {
                tiler.add(p).map_err(|e| Error::Other(format!("Invalid polygon in multipolygon: {}", e)))?;
            }
        }
        Geometry::GeometryCollection(gc) => {
            for geom in gc {
                match geom {
                    Geometry::Polygon(_) | Geometry::MultiPolygon(_) => {
                        feed_geometry(tiler, geom)?;
                    }
                    other => {
                        return Err(Error::Other(format!(
                            "Unsupported geometry type inside GeometryCollection: '{}'. \
                             Only Polygon and MultiPolygon are supported.",
                            geometry_type_name(&other)
                        )));
                    }
                }
            }
        }
        other => {
            return Err(Error::Other(format!(
                "Unsupported top-level geometry type: '{}'. \
                 Must be Polygon, MultiPolygon, or GeometryCollection.",
                geometry_type_name(&other)
            )));
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Base cell shuffle
// ---------------------------------------------------------------------------

/// Return a shuffled mapping of base cell indices 0–121.
/// This forms the most significant level of the GRTS address, ensuring that
/// sample coverage is balanced across the entire Earth surface.
fn shuffled_base_cells(seed_u64: u64) -> [u64; H3_NUM_BASE_CELLS] {
    let mut cells: [u64; H3_NUM_BASE_CELLS] = std::array::from_fn(|i| i as u64);
    let mut rng = Pcg64Mcg::seed_from_u64(seed_u64);
    cells.shuffle(&mut rng);
    cells
}


/// Return a human-readable name for a `geo::Geometry` variant.
fn geometry_type_name<T: geo::CoordNum>(geom: &Geometry<T>) -> &'static str {
    match geom {
        Geometry::Point(_)              => "Point",
        Geometry::Line(_)               => "Line",
        Geometry::LineString(_)         => "LineString",
        Geometry::Polygon(_)            => "Polygon",
        Geometry::MultiPoint(_)         => "MultiPoint",
        Geometry::MultiLineString(_)    => "MultiLineString",
        Geometry::MultiPolygon(_)       => "MultiPolygon",
        Geometry::GeometryCollection(_) => "GeometryCollection",
        Geometry::Rect(_)               => "Rect",
        Geometry::Triangle(_)           => "Triangle",
    }
}

// ---------------------------------------------------------------------------
// extendr module registration
// ---------------------------------------------------------------------------

extendr_module! {
    mod h3sampling;
    fn generate_grts_sample;
}
