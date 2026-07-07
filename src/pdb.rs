//! Fixed-column PDB ATOM parser producing the same atom universe biopython's
//! PDBParser hands to NeighborSearch: first model, highest-occupancy alternate
//! per disordered atom, standard amino-acid residues, hydrogens dropped.

use std::collections::HashMap;

use rsomics_common::Result;
use rsomics_common::error::RsomicsError;

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

const STANDARD_AA: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
];

fn is_standard_aa(resname: &str) -> bool {
    STANDARD_AA.contains(&resname)
}

/// biopython's PDBParser stores coordinates as 32-bit floats and only widens
/// them to f64 inside NeighborSearch's KDTree. Matching that f32 rounding is
/// what makes the inclusive `distance <= cutoff` boundary bit-exact against the
/// oracle, so coordinates are parsed through f32 here.
fn coord(field: &str) -> Result<f64> {
    let v: f32 = field.trim().parse()?;
    Ok(f64::from(v))
}

/// Element column 76-77 first, then the atom-name column fallback — mirrors how
/// biopython tags hydrogens (Atom.element, falling back to the name's leading
/// non-digit character).
fn is_hydrogen(line: &[u8], atom_name: &str) -> bool {
    let element = if line.len() >= 78 {
        std::str::from_utf8(&line[76..78])
            .unwrap_or("")
            .trim()
            .to_string()
    } else {
        String::new()
    };
    if !element.is_empty() {
        return element == "H" || element == "D";
    }
    let first = atom_name.chars().find(|c| !c.is_ascii_digit());
    matches!(first, Some('H') | Some('D'))
}

/// Occupancy (columns 55-60) picks the representative of a disordered atom.
/// biopython stores it as a plain float, so it is compared in f64 here — only
/// the strict `>` at tie boundaries is sensitive to the precision.
fn occupancy(line: &str, len: usize) -> Result<f64> {
    if len < 60 {
        return Ok(1.0);
    }
    let field = line[54..60].trim();
    if field.is_empty() {
        return Ok(1.0);
    }
    Ok(field.parse()?)
}

struct Selected {
    atom: Atom,
    occupancy: f64,
}

/// Parse the first model into the flat atom list NeighborSearch operates on.
///
/// Alternate locations collapse the way biopython's `DisorderedAtom` does:
/// atoms are grouped by (chain, residue, atom name) and the alternate with the
/// highest occupancy represents the group, regardless of its altloc label. A
/// strict occupancy tie keeps the first-encountered alternate, matching
/// biopython's `disordered_select` default. Each group holds the file position
/// of its first alternate, so the flat list order is unchanged for structures
/// without altlocs.
pub fn parse(text: &str) -> Result<Vec<Atom>> {
    let mut groups: Vec<Selected> = Vec::new();
    let mut index: HashMap<(char, i32, char, String), usize> = HashMap::new();
    let mut seen_model = false;

    for line in text.lines() {
        let bytes = line.as_bytes();
        if line.starts_with("ENDMDL") {
            break;
        }
        if line.starts_with("MODEL") {
            if seen_model {
                break;
            }
            seen_model = true;
            continue;
        }
        if !line.starts_with("ATOM") {
            continue;
        }
        if bytes.len() < 54 {
            return Err(RsomicsError::InvalidInput(format!(
                "truncated ATOM record: {line}"
            )));
        }

        let resname = line[17..20].trim().to_string();
        if !is_standard_aa(&resname) {
            continue;
        }
        let atom_name = line[12..16].trim().to_string();
        if is_hydrogen(bytes, &atom_name) {
            continue;
        }

        let chain = bytes[21] as char;
        let resseq: i32 = line[22..26].trim().parse()?;
        let icode = bytes[26] as char;

        let serial: i64 = line[6..11].trim().parse()?;
        let x = coord(&line[30..38])?;
        let y = coord(&line[38..46])?;
        let z = coord(&line[46..54])?;
        let occ = occupancy(line, bytes.len())?;

        let key = (chain, resseq, icode, atom_name.clone());
        let atom = Atom {
            chain,
            resseq,
            icode,
            resname,
            atom_name,
            serial,
            x,
            y,
            z,
        };
        match index.get(&key) {
            None => {
                index.insert(key, groups.len());
                groups.push(Selected {
                    atom,
                    occupancy: occ,
                });
            }
            Some(&gi) if occ > groups[gi].occupancy => {
                groups[gi] = Selected {
                    atom,
                    occupancy: occ,
                };
            }
            Some(_) => {}
        }
    }

    Ok(groups.into_iter().map(|g| g.atom).collect())
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
