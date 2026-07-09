//! Atom universe biopython's PDBParser hands to NeighborSearch, sourced from
//! the shared `rsomics-pdb-core` parser: first model, highest-occupancy
//! alternate per disordered atom, standard amino-acid residues, hydrogens
//! dropped.

use rsomics_common::Result;
use rsomics_pdb_core::{AltLocPolicy, ParseOptions, ResidueFilter, select_altloc};

pub struct Atom {
    pub chain: char,
    pub resseq: i32,
    pub icode: char,
    pub resname: String,
    pub atom_name: String,
    pub serial: i64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Parse the first model into the flat atom list NeighborSearch operates on:
/// standard amino acids only, hydrogens and HETATM dropped, and the
/// highest-occupancy alternate representing each disordered atom (ties keep the
/// first-seen alternate).
///
/// biopython's PDBParser stores coordinates as 32-bit floats and only widens
/// them to f64 inside NeighborSearch's KDTree. Matching that f32 rounding is
/// what makes the inclusive `distance <= cutoff` boundary bit-exact against the
/// oracle, so coordinates are rounded through f32 here.
pub fn parse(text: &str) -> Result<Vec<Atom>> {
    let opts = ParseOptions {
        include_hetatm: false,
        include_hydrogen: false,
        residue_filter: ResidueFilter::Standard20,
    };
    let atoms = select_altloc(
        rsomics_pdb_core::parse(text, &opts)?,
        AltLocPolicy::HighestOccupancy,
    );

    Ok(atoms
        .into_iter()
        .map(|a| Atom {
            chain: a.chain,
            resseq: a.resseq,
            icode: a.icode,
            resname: a.resname,
            atom_name: a.name,
            serial: a.serial,
            x: f64::from(a.coord[0] as f32),
            y: f64::from(a.coord[1] as f32),
            z: f64::from(a.coord[2] as f32),
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C
ATOM      3  H   THR A   1      16.000  12.000   4.000  1.00  0.00           H
HETATM    4  O   HOH A 100       1.000   2.000   3.000  1.00  0.00           O
ATOM      5  CB AALA B   2       1.000   2.000   3.000  0.50  0.00           C
ATOM      6  CB BALA B   2       9.000   9.000   9.000  0.50  0.00           C
";

    #[test]
    fn skips_hydrogen_hetatm_and_keeps_standard_ca() {
        let atoms = parse(SAMPLE).unwrap();
        let names: Vec<_> = atoms.iter().map(|a| a.atom_name.as_str()).collect();
        assert_eq!(
            names,
            vec!["N", "CA", "CB"],
            "H, HOH dropped; disordered CB collapsed to one"
        );
        assert_eq!(atoms[2].chain, 'B');
        // occupancy tie (0.50 == 0.50) keeps the first-seen alternate 'A'
        assert_eq!(atoms[2].x, 1.0, "tie keeps first-seen alternate");
    }

    #[test]
    fn coordinates_round_through_f32() {
        let atoms = parse(SAMPLE).unwrap();
        // 16.967 has no exact f64 form; biopython sees its f32 rounding
        assert_eq!(atoms[1].x, f64::from(16.967_f32));
    }

    #[test]
    fn highest_occupancy_alternate_wins() {
        // altloc A (occ 0.40) seen first, altloc B (occ 0.60) — biopython selects B
        let text = "\
ATOM      1  CA AALA A   1       0.000   0.000   0.000  0.40  0.00           C
ATOM      2  CA BALA A   1      10.000   0.000   0.000  0.60  0.00           C
";
        let atoms = parse(text).unwrap();
        assert_eq!(atoms.len(), 1, "the two alternates collapse to one");
        assert_eq!(atoms[0].x, 10.0, "highest-occupancy alternate B wins");
        assert_eq!(atoms[0].serial, 2, "selected serial is the B alternate's");
    }

    #[test]
    fn non_a_only_alternate_is_kept() {
        // an atom whose sole label is a non-'A' alternate must not be dropped
        let text = "\
ATOM      1  CA BALA A   1       5.000   0.000   0.000  1.00  0.00           C
ATOM      2  CB cALA A   1       6.000   0.000   0.000  1.00  0.00           C
ATOM      3  CG 2ALA A   1       7.000   0.000   0.000  1.00  0.00           C
";
        let atoms = parse(text).unwrap();
        let names: Vec<_> = atoms.iter().map(|a| a.atom_name.as_str()).collect();
        assert_eq!(
            names,
            vec!["CA", "CB", "CG"],
            "B-only, lowercase, and numeric alternates are all kept"
        );
    }
}
