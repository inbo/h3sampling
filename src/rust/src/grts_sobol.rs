// grts_sobol.rs — Exact-Area 1D-to-2D GRTS Routing
//
// -- The Algorithm
// This implements a GRTS spatial sampler without discrete grid approximations
// or rejection sampling.
//
// 1. 2D Sequence: Generates X and Y from independent Sobol dimensions (Van der
//    Corput and Joe-Kuo), applying independent Owen scrambles.
// 2. Morton Interleave: Interleaves the 32-bit X and Y into a 64-bit Z-curve
//    code, explicitly mapping the base-4 Owen-scrambled Halton sequence.
// 3. Area Targets: Converts the Morton code to a [0, 1) fraction and multiplies
//    by the exact total spherical area of the polygon.
// 4. Lazy Routing: Traverses a top-down quadtree. At each node, clips the
//    polygon geometry against the 4 child quadrants using Sutherland-Hodgman.
// 5. CDF Assignment: Partitions the sorted 1D points into the quadrants based
//    on the exact computed area CDF. Branches with 0 points are skipped entirely.
// 6. Prefix Stability: Re-sorts the final accepted points back to their original
//    sequence index to guarantee the Robertson master sample property.

use extendr_api::prelude::*;
use geo::algorithm::bounding_rect::BoundingRect;
use geo::{Coord, Geometry, LineString, Rect};
use geo_traits::to_geo::ToGeoGeometry;
use wkb::reader::Wkb;

const LON_MIN: f64 = -180.0;
const LON_MAX: f64 = 180.0;
const S_MIN: f64 = -1.0;
const S_MAX: f64 = 1.0;

const LEAF_DEPTH: u32 = 28;
const DENSIFY_DEG: f64 = 0.01;

// ---------------------------------------------------------------------------
// Coordinate transforms
// ---------------------------------------------------------------------------

#[inline]
fn lat_to_s(lat_deg: f64) -> f64 {
    lat_deg.to_radians().sin()
}

#[inline]
fn s_to_lat(s: f64) -> f64 {
    s.clamp(-1.0, 1.0).asin().to_degrees()
}

// ---------------------------------------------------------------------------
// Target Point Tracking
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct TargetPoint {
    id: usize,        // generation index (0..n-1); restored after routing for prefix stability
    grts_rank: usize, // position in area_target sort order = spatial GRTS sequence order
    area_target: f64,
    lon: f64,
    lat: f64,
}

// ---------------------------------------------------------------------------
// Exact Geometry & Clipping Structures
// ---------------------------------------------------------------------------

#[derive(Clone)]
struct SignedRing {
    coords: Vec<Coord>,
    sign: f64, // +1.0 for exteriors, -1.0 for holes
}

/// Densify a line string and map to (lon, sin(lat))
fn densify_and_transform_ring(ring: &LineString) -> Vec<Coord> {
    let coords = ring.0.as_slice();
    if coords.len() < 2 {
        return coords
            .iter()
            .map(|c| Coord {
                x: c.x,
                y: lat_to_s(c.y),
            })
            .collect();
    }

    let mut out: Vec<Coord> = Vec::with_capacity(coords.len() * 4);

    for window in coords.windows(2) {
        let (p0, p1) = (window[0], window[1]);
        let max_span = (p1.x - p0.x).abs().max((p1.y - p0.y).abs());

        let n_segs = if max_span > DENSIFY_DEG {
            (max_span / DENSIFY_DEG).ceil() as usize
        } else {
            1
        };

        for k in 0..n_segs {
            let t = k as f64 / n_segs as f64;
            let lon = p0.x + t * (p1.x - p0.x);
            let lat = p0.y + t * (p1.y - p0.y);
            out.push(Coord {
                x: lon,
                y: lat_to_s(lat),
            });
        }
    }

    let last = coords.last().unwrap();
    out.push(Coord {
        x: last.x,
        y: lat_to_s(last.y),
    });
    out
}

fn extract_signed_rings(geometry: &Geometry) -> extendr_api::Result<Vec<SignedRing>> {
    let mut rings = Vec::new();
    collect_rings(geometry, &mut rings)?;
    Ok(rings)
}

/// Recursively collect signed rings from any geometry type.
/// GeometryCollections are traversed recursively, so nested collections
/// (e.g. a GeometryCollection containing a MultiPolygon) are handled correctly.
fn collect_rings(geometry: &Geometry, rings: &mut Vec<SignedRing>) -> extendr_api::Result<()> {
    match geometry {
        Geometry::Polygon(p) => {
            rings.push(SignedRing {
                coords: densify_and_transform_ring(p.exterior()),
                sign: 1.0,
            });
            for hole in p.interiors() {
                rings.push(SignedRing {
                    coords: densify_and_transform_ring(hole),
                    sign: -1.0,
                });
            }
        }
        Geometry::MultiPolygon(mp) => {
            for p in &mp.0 {
                collect_rings(&Geometry::Polygon(p.clone()), rings)?;
            }
        }
        Geometry::GeometryCollection(gc) => {
            for g in gc.iter() {
                collect_rings(g, rings)?;
            }
        }
        _ => return Err(Error::Other("Must be a POLYGON or MULTIPOLYGON".into())),
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Sutherland-Hodgman Polygon Clipping & Area
// ---------------------------------------------------------------------------

#[inline]
fn inside(p: &Coord, edge: u8, val: f64) -> bool {
    match edge {
        0 => p.x >= val, // Left
        1 => p.x <= val, // Right
        2 => p.y >= val, // Bottom
        3 => p.y <= val, // Top
        _ => unreachable!(),
    }
}

#[inline]
fn intersect(p1: &Coord, p2: &Coord, edge: u8, val: f64) -> Coord {
    match edge {
        0 | 1 => Coord {
            x: val,
            y: p1.y + (val - p1.x) / (p2.x - p1.x) * (p2.y - p1.y),
        },
        2 | 3 => Coord {
            x: p1.x + (val - p1.y) / (p2.y - p1.y) * (p2.x - p1.x),
            y: val,
        },
        _ => unreachable!(),
    }
}

fn clip_edge(input: &[Coord], edge: u8, val: f64) -> Vec<Coord> {
    let mut out = Vec::with_capacity(input.len());
    if input.is_empty() {
        return out;
    }

    let mut prev = *input.last().unwrap();
    let mut prev_in = inside(&prev, edge, val);

    for &curr in input {
        let curr_in = inside(&curr, edge, val);
        if curr_in != prev_in {
            out.push(intersect(&prev, &curr, edge, val));
        }
        if curr_in {
            out.push(curr);
        }
        prev = curr;
        prev_in = curr_in;
    }
    out
}

fn clip_sh(input: &[Coord], xmin: f64, xmax: f64, ymin: f64, ymax: f64) -> Vec<Coord> {
    let p = clip_edge(input, 0, xmin);
    let p = clip_edge(&p, 1, xmax);
    let p = clip_edge(&p, 2, ymin);
    clip_edge(&p, 3, ymax)
}

fn shoelace_abs(coords: &[Coord]) -> f64 {
    if coords.len() < 3 {
        return 0.0;
    }
    let mut area = 0.0;
    for i in 0..coords.len() {
        let p1 = coords[i];
        let p2 = coords[(i + 1) % coords.len()];
        area += p1.x * p2.y - p2.x * p1.y;
    }
    (area / 2.0).abs()
}

fn calculate_clipped_area(rings: &[SignedRing], xmin: f64, xmax: f64, ymin: f64, ymax: f64) -> f64 {
    let mut total_area = 0.0;
    for ring in rings {
        let clipped = clip_sh(&ring.coords, xmin, xmax, ymin, ymax);
        if clipped.len() >= 3 {
            total_area += shoelace_abs(&clipped) * ring.sign;
        }
    }
    total_area.max(0.0) // Clamp to 0.0 to handle float anomalies
}

// ---------------------------------------------------------------------------
// 2D Sobol -> 1D Morton Sequence
// ---------------------------------------------------------------------------

#[inline]
fn lk_hash(mut x: u32, seed: u32) -> u32 {
    x ^= x.wrapping_mul(0x3d20_adea);
    x = x.wrapping_add(seed);
    x = x.wrapping_mul((seed >> 16) | 1);
    x ^= x.wrapping_mul(0x0552_6c56);
    x ^= x.wrapping_mul(0x53a2_2864);
    x
}

#[inline]
fn dim_seed(global_seed: u32, dim: usize) -> u32 {
    lk_hash(global_seed, dim as u32).wrapping_add(0x9e37_79b9)
}

#[inline]
fn owen_scramble(x: u32, seed: u32) -> u32 {
    let x = x.reverse_bits();
    let x = lk_hash(x, seed);
    x.reverse_bits()
}

const DIRECTION_VECTORS: [[u32; 32]; 2] = [
    // dim 0: Van der Corput
    [
        0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x08000000, 0x04000000, 0x02000000,
        0x01000000, 0x00800000, 0x00400000, 0x00200000, 0x00100000, 0x00080000, 0x00040000,
        0x00020000, 0x00010000, 0x00008000, 0x00004000, 0x00002000, 0x00001000, 0x00000800,
        0x00000400, 0x00000200, 0x00000100, 0x00000080, 0x00000040, 0x00000020, 0x00000010,
        0x00000008, 0x00000004, 0x00000002, 0x00000001,
    ],
    // dim 1: Joe-Kuo new-6
    [
        0x80000000, 0xc0000000, 0xa0000000, 0xf0000000, 0x88000000, 0xcc000000, 0xaa000000,
        0xff000000, 0x80800000, 0xc0c00000, 0xa0a00000, 0xf0f00000, 0x88880000, 0xcccc0000,
        0xaaaa0000, 0xffff0000, 0x80008000, 0xc000c000, 0xa000a000, 0xf000f000, 0x88008800,
        0xcc00cc00, 0xaa00aa00, 0xff00ff00, 0x80808080, 0xc0c0c0c0, 0xa0a0a0a0, 0xf0f0f0f0,
        0x88888888, 0xcccccccc, 0xaaaaaaaa, 0xffffffff,
    ],
];

#[inline]
fn sobol_raw(mut index: u32, dimension: usize) -> u32 {
    let v = &DIRECTION_VECTORS[dimension];
    let mut result = 0u32;
    let mut i = 0usize;
    while index > 0 {
        if index & 1 != 0 {
            result ^= v[i];
        }
        index >>= 1;
        i += 1;
    }
    result
}

/// Interleave the bits of X and Y into a 64-bit Morton code.
#[inline]
fn morton_interleave(x: u32, y: u32) -> u64 {
    let mut res = 0u64;
    for i in 0..32 {
        let shift = 31 - i;
        let bit_x = ((x >> shift) & 1) as u64;
        let bit_y = ((y >> shift) & 1) as u64;

        // Quadtree quadrant order: Y is the high bit, X is the low bit.
        // e.g., Quadrant 3 (Top-Right) -> X=1, Y=1 -> 0b11
        let out_shift = 62 - 2 * i;
        res |= bit_x << out_shift;
        res |= bit_y << (out_shift + 1);
    }
    res
}

// ---------------------------------------------------------------------------
// Quadtree Routing
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct CellBounds {
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
}

fn route(
    pts: &mut [TargetPoint],
    rings: &[SignedRing],
    cell: CellBounds,
    depth: u32,
    base_area: f64,
) {
    if pts.is_empty() {
        return;
    }

    if depth == LEAF_DEPTH {
        // Place all points at the leaf cell centre. At depth 28 the cell is
        // ~15 cm × ~5 cm, so the maximum placement error is ~8 cm — negligible
        // for any ecological application. Deterministic centre placement avoids
        // the within-cell jitter aliasing that arises when using additional Sobol
        // dimensions for sub-cell offsets (two indices sharing the same cell but
        // different low-order bits can produce colliding jitter values).
        // Two distinct area_targets cannot land in the same leaf cell because the
        // minimum Sobol inter-point gap (~2.3e-10) exceeds the leaf area fraction
        // (~1.4e-17), so pts always has exactly one element here in practice.
        let cx = (cell.xmin + cell.xmax) * 0.5;
        let cy = (cell.ymin + cell.ymax) * 0.5;
        for p in pts.iter_mut() {
            p.lon = cx;
            p.lat = cy;
        }
        return;
    }

    // Isolate rings that overlap this specific node
    let mut local_rings = Vec::new();
    for ring in rings {
        let clipped = clip_sh(&ring.coords, cell.xmin, cell.xmax, cell.ymin, cell.ymax);
        if clipped.len() >= 3 {
            local_rings.push(SignedRing {
                coords: clipped,
                sign: ring.sign,
            });
        }
    }

    if local_rings.is_empty() {
        return;
    }

    let xmid = (cell.xmin + cell.xmax) * 0.5;
    let ymid = (cell.ymin + cell.ymax) * 0.5;

    // Calculate exact area in the first 3 child quadrants.
    // We skip a3 because any remaining points automatically fall into quadrant 3!
    let a0 = calculate_clipped_area(&local_rings, cell.xmin, xmid, cell.ymin, ymid);
    let a1 = calculate_clipped_area(&local_rings, xmid, cell.xmax, cell.ymin, ymid);
    let a2 = calculate_clipped_area(&local_rings, cell.xmin, xmid, ymid, cell.ymax);

    // Absolute Area Thresholds (Z-curve order)
    let t0 = base_area + a0;
    let t1 = t0 + a1;
    let t2 = t1 + a2;

    // Partition points into the 4 child quadrants based on their area_target
    let mut i0 = 0;
    while i0 < pts.len() && pts[i0].area_target < t0 {
        i0 += 1;
    }
    let (pts0, rest) = pts.split_at_mut(i0);

    let mut i1 = 0;
    while i1 < rest.len() && rest[i1].area_target < t1 {
        i1 += 1;
    }
    let (pts1, rest) = rest.split_at_mut(i1);

    let mut i2 = 0;
    while i2 < rest.len() && rest[i2].area_target < t2 {
        i2 += 1;
    }
    let (pts2, pts3) = rest.split_at_mut(i2);

    // Route down branches using struct update syntax for cleaner bounds
    route(
        pts0,
        &local_rings,
        CellBounds {
            xmax: xmid,
            ymax: ymid,
            ..cell
        },
        depth + 1,
        base_area,
    );
    route(
        pts1,
        &local_rings,
        CellBounds {
            xmin: xmid,
            ymax: ymid,
            ..cell
        },
        depth + 1,
        t0,
    );
    route(
        pts2,
        &local_rings,
        CellBounds {
            xmax: xmid,
            ymin: ymid,
            ..cell
        },
        depth + 1,
        t1,
    );
    route(
        pts3,
        &local_rings,
        CellBounds {
            xmin: xmid,
            ymin: ymid,
            ..cell
        },
        depth + 1,
        t2,
    );
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

#[extendr]
pub fn grts_sobol_sample_from_wkb(
    wkb_bytes: &[u8],
    n: i32,
    seed: f64,
) -> extendr_api::Result<Robj> {
    let target_n = n as usize;
    if target_n == 0 {
        return Err(Error::Other("n must be >= 1".into()));
    }

    // Fold all 64 IEEE 754 bits of the f64 seed into a u32.
    // XOR-folding is preferred over truncation (seed as u64 as u32) because it
    // incorporates the exponent bits, ensuring nearby f64 values (e.g. 1.0 and 2.0)
    // produce different u32 seeds even when their low 32 bits are both zero.
    // seed = 0.0 gives seed_u32 = 0; dim_seed(0, d) = lk_hash(0, d) + golden_ratio,
    // which is non-zero for all d — seed=0 is a valid, non-degenerate choice.
    let seed_u32 = (seed.to_bits() ^ (seed.to_bits() >> 32)) as u32;
    let seed_x = dim_seed(seed_u32, 0);
    let seed_y = dim_seed(seed_u32, 1);

    // Decode WKB
    let wkb_obj =
        Wkb::try_new(wkb_bytes).map_err(|e| Error::Other(format!("WKB parse failed: {e}")))?;
    let geom_lonlat = wkb_obj.to_geometry();

    let bbox_lonlat = bounding_rect_geom(&geom_lonlat)
        .ok_or_else(|| Error::Other("Could not compute bounding box".into()))?;

    // Pre-calculate exact rings and total polygon area
    let signed_rings = extract_signed_rings(&geom_lonlat)?;
    let total_area = calculate_clipped_area(&signed_rings, LON_MIN, LON_MAX, S_MIN, S_MAX);

    if total_area <= 0.0 {
        return Err(Error::Other("Geometry has zero computable area".into()));
    }

    // 1D GRTS Generation via Interleaved 2D Sobol
    let mut targets = Vec::with_capacity(target_n);
    for i in 0..target_n {
        let raw_x = sobol_raw(i as u32, 0);
        let raw_y = sobol_raw(i as u32, 1);

        let sx = owen_scramble(raw_x, seed_x);
        let sy = owen_scramble(raw_y, seed_y);

        let m64 = morton_interleave(sx, sy);

        // Convert top 53 bits to an exact f64 fraction in [0, 1)
        let u = (m64 >> 11) as f64 / (1u64 << 53) as f64;

        targets.push(TargetPoint {
            id: i,
            grts_rank: 0, // filled after sorting below
            area_target: u * total_area,
            lon: 0.0,
            lat: 0.0,
        });
    }

    // Sort to enable fast routing partition, then record spatial rank.
    // grts_rank is the position in this sorted order — it defines the GRTS
    // sequence: rank 0 is the first point visited spatially, rank 1 the second,
    // etc. This ordering is useful for sequential field sampling (visit sites
    // in rank order to maintain spatial balance as sampling is extended or stopped).
    targets.sort_by(|a, b| a.area_target.partial_cmp(&b.area_target).unwrap());
    for (rank, t) in targets.iter_mut().enumerate() {
        t.grts_rank = rank;
    }

    // Route points through the quadtree
    let root_cell = CellBounds {
        xmin: LON_MIN,
        xmax: LON_MAX,
        ymin: S_MIN,
        ymax: S_MAX,
    };
    route(&mut targets, &signed_rings, root_cell, 1, 0.0);

    // Un-sort points back to their exact sequence generation order (prefix stability).
    targets.sort_by_key(|t| t.id);

    let lons: Vec<f64> = targets.iter().map(|t| t.lon).collect();
    let lats: Vec<f64> = targets.iter().map(|t| s_to_lat(t.lat)).collect();
    // grts_rank gives the spatial ordering: visit sites in ascending rank order
    // to maintain spatial balance when the sample is extended or truncated in the field.
    let grts_ranks: Vec<i32> = targets.iter().map(|t| t.grts_rank as i32).collect();

    let df = data_frame!(lon = lons, lat = lats, grts_rank = grts_ranks);

    let mut df_robj: Robj = df;
    df_robj.set_attrib("seed", seed)?;
    df_robj.set_attrib("n", n)?;
    df_robj.set_attrib(
        "bbox",
        r!([
            bbox_lonlat.min().x,
            bbox_lonlat.min().y,
            bbox_lonlat.max().x,
            bbox_lonlat.max().y
        ]),
    )?;
    df_robj.set_attrib("sobol_scanned", target_n as f64)?;
    df_robj.set_attrib("acceptance_rate", 1.0)?; // Zero rejection algorithm

    Ok(df_robj)
}

fn bounding_rect_geom(geometry: &Geometry) -> Option<Rect> {
    match geometry {
        Geometry::Polygon(p) => p.bounding_rect(),
        Geometry::MultiPolygon(m) => m.bounding_rect(),
        Geometry::GeometryCollection(gc) => {
            gc.iter().filter_map(bounding_rect_geom).reduce(|a, b| {
                Rect::new(
                    Coord {
                        x: a.min().x.min(b.min().x),
                        y: a.min().y.min(b.min().y),
                    },
                    Coord {
                        x: a.max().x.max(b.max().x),
                        y: a.max().y.max(b.max().y),
                    },
                )
            })
        }
        _ => None,
    }
}

extendr_module! {
    mod grts_sobol;
    fn grts_sobol_sample_from_wkb;
}
