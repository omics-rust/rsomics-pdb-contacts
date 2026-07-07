//! Contact-pair extraction matching Bio.PDB.NeighborSearch.search_all.
//!
//! search_all(cutoff, level='R') returns the residue pairs that have at least
//! one atom pair within `cutoff`; level='A' returns the atom pairs themselves.
//! The boundary is inclusive (distance <= cutoff). Self-pairs never occur and,
//! at residue level, atom pairs inside one residue are dropped. We additionally
//! expose --ca-only (restrict the atom universe to Cα) and --min-seq-sep (skip
//! intra-chain pairs whose residue-sequence indices are too close) as
//! conventions layered on top of the same predicate.

use rayon::prelude::*;
use rsomics_common::Result;
use serde::Serialize;

use crate::cell_list::{Point, neighbor_pairs};
use crate::pdb::Atom;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Level {
    Residue,
    Atom,
}

/// Identity tuple ordered exactly like biopython's residue full_id[1:] for a
/// single model: (chain, hetfield, resseq, icode). ATOM records carry a blank
/// hetfield, so it sorts ahead of any HETATM residue uniformly.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
struct ResId {
    chain: char,
    resseq: i32,
    icode: char,
}

impl ResId {
    fn of(a: &Atom) -> Self {
        Self {
            chain: a.chain,
            resseq: a.resseq,
            icode: a.icode,
        }
    }
}

#[derive(Clone, Serialize)]
pub struct ResiduePair {
    pub chain_a: char,
    pub resseq_a: i32,
    pub icode_a: char,
    pub resname_a: String,
    pub chain_b: char,
    pub resseq_b: i32,
    pub icode_b: char,
    pub resname_b: String,
}

#[derive(Clone, Serialize)]
pub struct AtomPair {
    pub chain_a: char,
    pub resseq_a: i32,
    pub icode_a: char,
    pub resname_a: String,
    pub atom_a: String,
    pub chain_b: char,
    pub resseq_b: i32,
    pub icode_b: char,
    pub resname_b: String,
    pub atom_b: String,
}

#[derive(Serialize)]
pub enum ContactSet {
    Residue(Vec<ResiduePair>),
    Atom(Vec<AtomPair>),
}

/// Same-chain pair too close in sequence to count, given `min_seq_sep`. Matches
/// biopython-side filtering only on the residue sequence index (icode ignored,
/// as the convention is defined on resseq).
fn seq_excluded(a: &Atom, b: &Atom, min_seq_sep: i32) -> bool {
    min_seq_sep > 0 && a.chain == b.chain && (a.resseq - b.resseq).abs() < min_seq_sep
}

fn select_atoms(atoms: &[Atom], ca_only: bool) -> Vec<usize> {
    atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| !ca_only || a.atom_name == "CA")
        .map(|(i, _)| i)
        .collect()
}

pub fn contacts(
    atoms: &[Atom],
    cutoff: f64,
    level: Level,
    ca_only: bool,
    min_seq_sep: i32,
) -> Result<ContactSet> {
    let selected = select_atoms(atoms, ca_only);
    let points: Vec<Point> = selected
        .iter()
        .map(|&i| Point {
            x: atoms[i].x,
            y: atoms[i].y,
            z: atoms[i].z,
        })
        .collect();

    let raw = neighbor_pairs(&points, cutoff)?;

    Ok(match level {
        Level::Atom => ContactSet::Atom(atom_pairs(atoms, &selected, &raw, min_seq_sep)),
        Level::Residue => ContactSet::Residue(residue_pairs(atoms, &selected, &raw, min_seq_sep)),
    })
}

fn atom_pairs(
    atoms: &[Atom],
    selected: &[usize],
    raw: &[(usize, usize)],
    min_seq_sep: i32,
) -> Vec<AtomPair> {
    let mut pairs: Vec<(ResId, i64, ResId, i64)> = raw
        .par_iter()
        .filter_map(|&(pi, pj)| {
            let (ia, ja) = (selected[pi], selected[pj]);
            let (a, b) = (&atoms[ia], &atoms[ja]);
            if seq_excluded(a, b, min_seq_sep) {
                return None;
            }
            let ka = (ResId::of(a), a.serial);
            let kb = (ResId::of(b), b.serial);
            let ((ra, sa), (rb, sb)) = if ka <= kb { (ka, kb) } else { (kb, ka) };
            Some((ra, sa, rb, sb))
        })
        .collect();
    pairs.par_sort_unstable();

    let by_serial: std::collections::HashMap<i64, usize> =
        selected.iter().map(|&i| (atoms[i].serial, i)).collect();

    pairs
        .into_iter()
        .map(|(_, sa, _, sb)| {
            let a = &atoms[by_serial[&sa]];
            let b = &atoms[by_serial[&sb]];
            AtomPair {
                chain_a: a.chain,
                resseq_a: a.resseq,
                icode_a: a.icode,
                resname_a: a.resname.clone(),
                atom_a: a.atom_name.clone(),
                chain_b: b.chain,
                resseq_b: b.resseq,
                icode_b: b.icode,
                resname_b: b.resname.clone(),
                atom_b: b.atom_name.clone(),
            }
        })
        .collect()
}

fn residue_pairs(
    atoms: &[Atom],
    selected: &[usize],
    raw: &[(usize, usize)],
    min_seq_sep: i32,
) -> Vec<ResiduePair> {
    let mut keyed: Vec<(ResId, ResId)> = raw
        .par_iter()
        .filter_map(|&(pi, pj)| {
            let (ia, ja) = (selected[pi], selected[pj]);
            let (a, b) = (&atoms[ia], &atoms[ja]);
            let ra = ResId::of(a);
            let rb = ResId::of(b);
            if ra == rb {
                return None;
            }
            if seq_excluded(a, b, min_seq_sep) {
                return None;
            }
            Some(if ra < rb { (ra, rb) } else { (rb, ra) })
        })
        .collect();
    keyed.par_sort_unstable();
    keyed.dedup();

    let resname: std::collections::HashMap<(char, i32, char), String> = atoms
        .iter()
        .map(|a| ((a.chain, a.resseq, a.icode), a.resname.clone()))
        .collect();
    let name = |r: &ResId| resname[&(r.chain, r.resseq, r.icode)].clone();

    keyed
        .into_iter()
        .map(|(ra, rb)| ResiduePair {
            chain_a: ra.chain,
            resseq_a: ra.resseq,
            icode_a: ra.icode,
            resname_a: name(&ra),
            chain_b: rb.chain,
            resseq_b: rb.resseq,
            icode_b: rb.icode,
            resname_b: name(&rb),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(chain: char, resseq: i32, name: &str, serial: i64, x: f64) -> Atom {
        Atom {
            chain,
            resseq,
            icode: ' ',
            resname: "ALA".to_string(),
            atom_name: name.to_string(),
            serial,
            x,
            y: 0.0,
            z: 0.0,
        }
    }

    fn res_pairs(set: ContactSet) -> Vec<ResiduePair> {
        match set {
            ContactSet::Residue(v) => v,
            ContactSet::Atom(_) => panic!("expected residue set"),
        }
    }

    fn atom_pairs_of(set: ContactSet) -> Vec<AtomPair> {
        match set {
            ContactSet::Atom(v) => v,
            ContactSet::Residue(_) => panic!("expected atom set"),
        }
    }

    #[test]
    fn intra_residue_atom_pairs_drop_at_residue_level() {
        // two atoms of one residue, 1.0 apart, plus a far residue
        let atoms = vec![
            atom('A', 1, "N", 1, 0.0),
            atom('A', 1, "CA", 2, 1.0),
            atom('A', 2, "N", 3, 2.0),
        ];
        let pairs = res_pairs(contacts(&atoms, 3.0, Level::Residue, false, 0).unwrap());
        // only residue 1 <-> residue 2; the intra-residue 1-1 pair is dropped
        assert_eq!(pairs.len(), 1);
        assert_eq!((pairs[0].resseq_a, pairs[0].resseq_b), (1, 2));
    }

    #[test]
    fn atom_level_keeps_intra_residue_pairs() {
        let atoms = vec![atom('A', 1, "N", 1, 0.0), atom('A', 1, "CA", 2, 1.0)];
        let pairs = atom_pairs_of(contacts(&atoms, 3.0, Level::Atom, false, 0).unwrap());
        assert_eq!(pairs.len(), 1);
        assert_eq!(
            (pairs[0].atom_a.as_str(), pairs[0].atom_b.as_str()),
            ("N", "CA")
        );
    }

    #[test]
    fn min_seq_sep_skips_close_intrachain_residues() {
        let atoms = vec![
            atom('A', 1, "CA", 1, 0.0),
            atom('A', 2, "CA", 2, 1.0),
            atom('A', 5, "CA", 3, 2.0),
        ];
        // sep<3 dropped: 1-2 (|1-2|=1) and 2-5? no |2-5|=3 ok ; 1-5 |4| ok
        let pairs = res_pairs(contacts(&atoms, 3.0, Level::Residue, false, 3).unwrap());
        let got: Vec<(i32, i32)> = pairs.iter().map(|p| (p.resseq_a, p.resseq_b)).collect();
        assert!(
            !got.contains(&(1, 2)),
            "1-2 within seq window must be skipped"
        );
        assert!(got.contains(&(2, 5)));
        assert!(got.contains(&(1, 5)));
    }

    #[test]
    fn ca_only_restricts_universe() {
        let atoms = vec![
            atom('A', 1, "CA", 1, 0.0),
            atom('A', 1, "CB", 2, 0.5),
            atom('A', 2, "CA", 3, 1.0),
        ];
        let pairs = atom_pairs_of(contacts(&atoms, 3.0, Level::Atom, true, 0).unwrap());
        // only the two CA atoms remain -> exactly one pair
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].atom_a, "CA");
        assert_eq!(pairs[0].atom_b, "CA");
    }

    #[test]
    fn pair_order_is_canonical_smaller_residue_first() {
        let atoms = vec![atom('B', 9, "CA", 2, 0.0), atom('A', 1, "CA", 1, 1.0)];
        let pairs = res_pairs(contacts(&atoms, 3.0, Level::Residue, false, 0).unwrap());
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].chain_a, 'A');
        assert_eq!(pairs[0].chain_b, 'B');
    }
}
