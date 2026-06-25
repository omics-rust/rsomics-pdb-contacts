//! Uniform-grid (cell-list) all-pairs neighbor search. Cells are `cutoff` wide,
//! so every pair within `cutoff` lies in adjacent cells and the 3x3x3 stencil
//! around each cell finds them. This is the same spatial-index idea freesasa's
//! SASA neighbor build uses, reimplemented here for the contact predicate.
//!
//! The contact test is `dx*dx + dy*dy + dz*dz <= cutoff*cutoff`, computed in
//! f64, matching biopython's kdtrees.c `KDTree_test_neighbors`
//! (`r <= self->_neighbor_radius_sq`) bit-for-bit on the boundary.

use rayon::prelude::*;

pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// All unordered index pairs `(i, j)` with `i < j` whose euclidean distance is
/// `<= cutoff`. Output order is unspecified; callers sort.
pub fn neighbor_pairs(points: &[Point], cutoff: f64) -> Vec<(usize, usize)> {
    let n = points.len();
    if n < 2 {
        return Vec::new();
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

    (0..nz)
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
        assert_eq!(sorted(neighbor_pairs(&pts, 2.0)), vec![(0, 1)]);
        assert_eq!(
            sorted(neighbor_pairs(&pts, 6.0)),
            vec![(0, 1), (0, 2), (1, 2)]
        );
    }

    #[test]
    fn boundary_is_inclusive() {
        let pts = vec![p(0.0, 0.0, 0.0), p(3.0, 4.0, 0.0)];
        assert_eq!(
            neighbor_pairs(&pts, 5.0).len(),
            1,
            "distance == cutoff included"
        );
        assert!(neighbor_pairs(&pts, 4.999_999).is_empty());
    }

    #[test]
    fn no_self_pairs_and_empty_below_two_points() {
        assert!(neighbor_pairs(&[], 5.0).is_empty());
        assert!(neighbor_pairs(&[p(0.0, 0.0, 0.0)], 5.0).is_empty());
    }
}
