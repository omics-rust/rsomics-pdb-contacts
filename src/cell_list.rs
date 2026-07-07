//! Uniform-grid (cell-list) all-pairs neighbor search. Cells are `cutoff` wide,
//! so every pair within `cutoff` lies in adjacent cells and the 3x3x3 stencil
//! around each cell finds them. This is the same spatial-index idea freesasa's
//! SASA neighbor build uses, reimplemented here for the contact predicate.
//!
//! The contact test is `dx*dx + dy*dy + dz*dz <= cutoff*cutoff`, computed in
//! f64, matching biopython's kdtrees.c `KDTree_test_neighbors`
//! (`r <= self->_neighbor_radius_sq`) bit-for-bit on the boundary.

use rayon::prelude::*;
use rsomics_common::Result;
use rsomics_common::error::RsomicsError;

pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// All unordered index pairs `(i, j)` with `i < j` whose euclidean distance is
/// `<= cutoff`. Output order is unspecified; callers sort.
///
/// A non-positive (or NaN) cutoff is rejected loudly, matching biopython's
/// `NeighborSearch.search_all` (`ValueError: Radius must be positive.`). A valid
/// but tiny cutoff would blow the uniform grid's cell count past any useful
/// bound; past a cap the grid is abandoned for the exact O(n^2) scan, which is
/// bounded for the moderate atom counts of a PDB and yields the same pair set.
pub fn neighbor_pairs(points: &[Point], cutoff: f64) -> Result<Vec<(usize, usize)>> {
    if cutoff.is_nan() || cutoff <= 0.0 {
        return Err(RsomicsError::InvalidInput(
            "radius must be positive".to_string(),
        ));
    }
    let n = points.len();
    if n < 2 {
        return Ok(Vec::new());
    }
    let cutoff_sq = cutoff * cutoff;
    let cell_size = cutoff;

    let (mut min_x, mut min_y, mut min_z) = (f64::MAX, f64::MAX, f64::MAX);
    let (mut max_x, mut max_y, mut max_z) = (f64::MIN, f64::MIN, f64::MIN);
    for p in points {
        min_x = min_x.min(p.x);
        min_y = min_y.min(p.y);
        min_z = min_z.min(p.z);
        max_x = max_x.max(p.x);
        max_y = max_y.max(p.y);
        max_z = max_z.max(p.z);
    }

    let dim = |lo: f64, hi: f64| -> i64 { (((hi - lo) / cell_size).floor() as i64) + 1 };
    let nx = dim(min_x, max_x);
    let ny = dim(min_y, max_y);
    let nz = dim(min_z, max_z);

    let total_cells = (nx as i128)
        .saturating_mul(ny as i128)
        .saturating_mul(nz as i128);
    let cap = (n as i128).saturating_mul(n as i128).max(1 << 20);
    if total_cells > cap {
        return Ok(brute_force(points, cutoff_sq));
    }

    let cell_of = |p: &Point| -> (i64, i64, i64) {
        (
            (((p.x - min_x) / cell_size).floor() as i64).clamp(0, nx - 1),
            (((p.y - min_y) / cell_size).floor() as i64).clamp(0, ny - 1),
            (((p.z - min_z) / cell_size).floor() as i64).clamp(0, nz - 1),
        )
    };
    let idx = |cx: i64, cy: i64, cz: i64| -> usize { ((cz * ny + cy) * nx + cx) as usize };

    let mut cells: Vec<Vec<usize>> = vec![Vec::new(); (nx * ny * nz) as usize];
    for (i, p) in points.iter().enumerate() {
        let (cx, cy, cz) = cell_of(p);
        cells[idx(cx, cy, cz)].push(i);
    }

    let within = |ia: usize, ja: usize| -> bool {
        let dx = points[ja].x - points[ia].x;
        let dy = points[ja].y - points[ia].y;
        let dz = points[ja].z - points[ia].z;
        dx * dx + dy * dy + dz * dz <= cutoff_sq
    };

    Ok((0..nz)
        .into_par_iter()
        .flat_map(|cz| {
            let mut local: Vec<(usize, usize)> = Vec::new();
            for cy in 0..ny {
                for cx in 0..nx {
                    let here = &cells[idx(cx, cy, cz)];
                    if here.is_empty() {
                        continue;
                    }
                    for dz in -1..=1 {
                        for dy in -1..=1 {
                            for dx in -1..=1 {
                                let (ox, oy, oz) = (cx + dx, cy + dy, cz + dz);
                                if ox < 0 || oy < 0 || oz < 0 || ox >= nx || oy >= ny || oz >= nz {
                                    continue;
                                }
                                for &ia in here {
                                    for &ja in &cells[idx(ox, oy, oz)] {
                                        if ia >= ja {
                                            continue;
                                        }
                                        if within(ia, ja) {
                                            local.push((ia, ja));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            local
        })
        .collect())
}

/// Exact all-pairs scan, the bounded fallback when a tiny cutoff makes the
/// uniform grid pathological. Same predicate and pair set as the grid path.
fn brute_force(points: &[Point], cutoff_sq: f64) -> Vec<(usize, usize)> {
    let n = points.len();
    (0..n)
        .into_par_iter()
        .flat_map(|i| {
            let mut local = Vec::new();
            for j in (i + 1)..n {
                let dx = points[j].x - points[i].x;
                let dy = points[j].y - points[i].y;
                let dz = points[j].z - points[i].z;
                if dx * dx + dy * dy + dz * dz <= cutoff_sq {
                    local.push((i, j));
                }
            }
            local
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64, z: f64) -> Point {
        Point { x, y, z }
    }

    fn sorted(mut v: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
        v.sort_unstable();
        v
    }

    #[test]
    fn finds_pairs_within_cutoff_only() {
        let pts = vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(5.0, 0.0, 0.0)];
        assert_eq!(sorted(neighbor_pairs(&pts, 2.0).unwrap()), vec![(0, 1)]);
        assert_eq!(
            sorted(neighbor_pairs(&pts, 6.0).unwrap()),
            vec![(0, 1), (0, 2), (1, 2)]
        );
    }

    #[test]
    fn boundary_is_inclusive() {
        let pts = vec![p(0.0, 0.0, 0.0), p(3.0, 4.0, 0.0)];
        assert_eq!(
            neighbor_pairs(&pts, 5.0).unwrap().len(),
            1,
            "distance == cutoff included"
        );
        assert!(neighbor_pairs(&pts, 4.999_999).unwrap().is_empty());
    }

    #[test]
    fn no_self_pairs_and_empty_below_two_points() {
        assert!(neighbor_pairs(&[], 5.0).unwrap().is_empty());
        assert!(neighbor_pairs(&[p(0.0, 0.0, 0.0)], 5.0).unwrap().is_empty());
    }

    #[test]
    fn rejects_nonpositive_cutoff() {
        let pts = vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0)];
        for bad in [0.0, -1.0, f64::NAN] {
            let err = neighbor_pairs(&pts, bad).unwrap_err();
            assert!(
                err.to_string().contains("radius must be positive"),
                "cutoff {bad} must be rejected loudly, got: {err}"
            );
        }
        // rejected even with too few points to form a pair
        assert!(neighbor_pairs(&[], 0.0).is_err());
    }

    #[test]
    fn tiny_cutoff_falls_back_without_overflow() {
        // a 1e6-wide span with a 0.01 cutoff would need ~1e8 grid cells; the
        // fallback must handle it and still find the coincident pair exactly.
        let pts = vec![p(0.0, 0.0, 0.0), p(0.0, 0.0, 0.0), p(1.0e6, 0.0, 0.0)];
        assert_eq!(sorted(neighbor_pairs(&pts, 0.01).unwrap()), vec![(0, 1)]);
    }
}
