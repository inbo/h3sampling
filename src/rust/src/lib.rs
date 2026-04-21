use extendr_api::prelude::*;
use h3o::{CellIndex, Resolution};
use h3o::geom::{TilerBuilder, ContainmentMode};
use geo::Geometry;
use geozero::wkb::Wkb;
use geozero::ToGeo;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
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
/// * `wkb_bytes`      – Study area geometry encoded as Well-Known Binary.
/// * `res`            – H3 resolution (0 = coarsest, 15 = finest).
/// * `containment`    – One of `"intersect"`, `"centroid"`, `"boundary"`, `"covers"`.
/// * `n`              – Desired sample size.
/// * `global_seed`    – Reproducibility seed (passed as f64 from R; values above
///                      2^53 lose precision — prefer seeds within that range).
/// * `area_correction`– If `false` (default), all H3 cells have equal inclusion
///                      probability pi_i = n / N. Cell areas vary by up to ~1.9x
///                      across the icosahedral projection; this is ignored in the
///                      sampling step but can be corrected at the estimation stage
///                      using the returned `area_m2` and `ip` columns.
///                      If `true`, rejection sampling is applied so that each cell
///                      is accepted with probability area_i / max_area, making the
///                      effective inclusion probability proportional to area:
///                      pi_i = n * area_i / sum(area). A pre-pass over the coverage
///                      is required to find the maximum cell area, so this mode
///                      incurs roughly 2x the runtime of the default mode.
///
/// # Returns
/// A `data.frame` with four columns ordered by ascending sort key (GRTS order):
/// * `cell`     – H3 cell identifier as a lowercase hexadecimal string.
/// * `sort_key` – GRTS sort key as a zero-padded 16-character lowercase hex
///                string. Zero-padding ensures correct lexicographic ordering
///                on the R side (e.g. for merging tile results).
/// * `area_m2` – True area of the sampled cell in m².
/// * `ip`       – Inclusion probability pi_i of the sampled cell:
///                  area_correction = false: pi_i = n / N  (constant)
///                  area_correction = true:  pi_i = n * area_i / sum(area)
///
/// Design-based estimators using the returned columns:
///   HT total:   T_hat   = sum(y_i * area_m2_i / ip_i)
///   Hajek mean: mu_Hajek = sum(y_i * area_m2_i / ip_i) / sum(area_m2_i / ip_i)
///
/// The returned `data.frame` also carries the following design attributes:
/// * `n_cells`        – Total H3 cells N covering the study area.
/// * `sum_area_m2`    – Sum of all cell areas in m².
/// * `n`              – Requested sample size.
/// * `resolution`     – H3 resolution used.
/// * `containment`    – Containment mode used.
/// * `area_correction`– Whether area-proportional sampling was applied.
#[extendr]
fn generate_grts_sample(
    wkb_bytes: &[u8],
    res: i32,
    containment: &str,
    n: i32,
    global_seed: f64,
    area_correction: bool,
) -> extendr_api::Result<Robj> {
    // NOTE: f64 -> u64 truncates; seeds above 2^53 lose precision.
    let seed_u64 = global_seed as u64;
    let target_n = n as usize;

    // 1. Decode WKB and resolve tiler settings ------------------------------
    let geometry = decode_wkb(wkb_bytes)?;
    let containment_mode = parse_containment(containment)?;
    let resolution = Resolution::try_from(res as u8)
        .map_err(|_| Error::Other(format!("Invalid H3 resolution: {}. Must be 0-15.", res)))?;

    // 2. Area-correction pre-pass ------------------------------------------
    // When area_correction is requested we need the maximum cell area at this
    // resolution before the main loop so that acceptance probabilities stay
    // in [0, 1]. The geometry is cloned to feed a second tiler.
    //
    // Mathematical note: the acceptance constant (max_area) cancels in the
    // inclusion probability formula. pi_i = n * area_i / sum(area) holds for
    // ANY constant C >= max_area used as the denominator in the acceptance
    // step. What matters is that C >= max_area to prevent acceptance
    // probabilities exceeding 1.
    let max_area_m2: f64 = if area_correction {
        let mut pre_tiler = TilerBuilder::new(resolution)
            .containment_mode(containment_mode)
            .build();
        feed_geometry(&mut pre_tiler, geometry.clone())?;
        pre_tiler
            .into_coverage()
            .map(|cell| cell.area_m2())
            .fold(0.0_f64, f64::max)
    } else {
        0.0 // unused in the equal-probability path
    };

    // 3. Build the main tiler ----------------------------------------------
    let mut tiler = TilerBuilder::new(resolution)
        .containment_mode(containment_mode)
        .build();
    feed_geometry(&mut tiler, geometry)?;

    // 4. Shuffle base cells (top level of the GRTS address) ----------------
    // This randomises which continental regions receive low sort keys,
    // ensuring global spatial balance across the entire Earth surface.
    let base_cells = shuffled_base_cells(seed_u64);

    // 5. Stream cells through the heap reservoir ----------------------------
    // The heap stores (sort_key, h3_index, area_bits) tuples. BinaryHeap is a
    // max-heap so peek() always returns the *largest* key — exactly what we
    // want to evict when a smaller key arrives. area_bits is the cell area
    // (m2) bit-cast to u64 so it can be stored without a separate Vec.
    let mut heap: BinaryHeap<(u64, u64, u64)> = BinaryHeap::with_capacity(target_n + 1);

    // Permutation cache: siblings share one RNG initialisation.
    let mut perm_cache: HashMap<u64, [u64; 7]> = HashMap::new();

    // Rejection RNG — seeded independently from the permutation RNG so that
    // the acceptance decision for one cell cannot affect the spatial ordering
    // of any other. The XOR constant decorrelates the two streams while still
    // being derived from the same user-facing seed.
    let mut rejection_rng = Pcg64Mcg::seed_from_u64(seed_u64 ^ 0xf0cacc1a);

    let mut total_cells: usize = 0;
    // Accumulate area over the FULL coverage regardless of acceptance, since
    // sum(area) — the denominator in pi_i — refers to the whole study area.
    let mut sum_all_areas: f64 = 0.0;

    for cell in tiler.into_coverage() {
        total_cells += 1;
        let cell_area = cell.area_m2();
        sum_all_areas += cell_area;

        // Rejection step (area_correction = true only) ---------------------
        // Accept cell with probability area_i / max_area. This thins the
        // coverage so that larger cells survive more often, making the
        // effective sampling proportional to physical area.
        if area_correction && rejection_rng.gen::<f64>() >= cell_area / max_area_m2 {
            continue;
        }

        let sort_key = compute_sort_key(cell, res, seed_u64, &base_cells, &mut perm_cache)?;
        reservoir_push(&mut heap, target_n, sort_key, u64::from(cell), cell_area);
    }

    // 6. Validate sample size -----------------------------------------------
    // Check heap size rather than total_cells: with area_correction the
    // eligible pool is smaller than the full coverage.
    if heap.len() < target_n {
        return Err(Error::Other(format!(
            "Target n ({}) exceeds the number of eligible cells ({}{}).              Try a higher resolution or a larger study area.",
            target_n,
            heap.len(),
            if area_correction {
                format!(" accepted from {} total after area correction", total_cells)
            } else {
                String::new()
            }
        )));
    }

    // 7. Extract results in GRTS order --------------------------------------
    // BinaryHeap::into_sorted_vec() drains the max-heap in ascending order,
    // which corresponds directly to the GRTS visiting sequence.
    let sorted: Vec<(u64, u64, u64)> = heap.into_sorted_vec();

    let cells: Vec<String> = sorted
        .iter()
        .map(|&(_, idx, _)| format!("{:x}", idx))
        .collect();

    let sort_keys: Vec<String> = sorted
        .iter()
        .map(|&(key, _, _)| format!("{:016x}", key))
        .collect();

    let areas_m2: Vec<f64> = sorted
        .iter()
        .map(|&(_, _, area_bits)| f64::from_bits(area_bits))
        .collect();

    // Inclusion probabilities ----------------------------------------------
    //   area_correction = false: pi_i = n / N  (equal for all cells)
    //   area_correction = true:  pi_i = n * area_i / sum(area)
    //
    // HT estimator in both cases: T_hat = sum(y_i * area_m2_i / ip_i)
    let n_f64 = target_n as f64;
    let ip: Vec<f64> = if area_correction {
        sorted
            .iter()
            .map(|&(_, _, area_bits)| {
                let area = f64::from_bits(area_bits);
                n_f64 * area / sum_all_areas
            })
            .collect()
    } else {
        vec![n_f64 / total_cells as f64; target_n]
    };

    // 8. Build the data.frame and attach design attributes ------------------
    let df = data_frame!(
        cell     = cells,
        sort_key = sort_keys,
        area_m2  = areas_m2,
        ip       = ip
    );

    let mut df_robj: Robj = df.into();
    df_robj.set_attrib("n_cells",         total_cells as i32)?;
    df_robj.set_attrib("sum_area_m2",    sum_all_areas)?;
    df_robj.set_attrib("n",               n)?;
    df_robj.set_attrib("resolution",      res)?;
    df_robj.set_attrib("containment",     containment)?;
    df_robj.set_attrib("area_correction", area_correction)?;

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
