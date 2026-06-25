//! Fixed-column PDB ATOM parser producing the same atom universe biopython's
//! PDBParser hands to NeighborSearch: first model, first altloc per residue,
//! standard amino-acid residues, hydrogens dropped.

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

/// Resolve altlocs the way biopython's default PDBParser does: the first altloc
/// character seen for a residue's atom wins; later differing altlocs of the same
/// atom are discarded.
struct AltlocFilter {
    last_res: Option<(char, i32, char)>,
    seen: Vec<(String, char)>,
}

impl AltlocFilter {
    fn new() -> Self {
        Self {
            last_res: None,
            seen: Vec::new(),
        }
    }

    fn accept(&mut self, key: (char, i32, char), atom_name: &str, altloc: char) -> bool {
        if self.last_res != Some(key) {
            self.last_res = Some(key);
            self.seen.clear();
        }
        if altloc == ' ' {
            return true;
        }
        if let Some((_, kept)) = self.seen.iter().find(|(n, _)| n == atom_name) {
            return *kept == altloc;
        }
        self.seen.push((atom_name.to_string(), altloc));
        true
    }
}

/// Parse the first model into the flat atom list NeighborSearch operates on.
pub fn parse(text: &str) -> Result<Vec<Atom>> {
    let mut atoms = Vec::new();
    let mut seen_model = false;
    let mut altloc = AltlocFilter::new();

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
        let alt = bytes[16] as char;
        if !altloc.accept((chain, resseq, icode), &atom_name, alt) {
            continue;
        }

        let serial: i64 = line[6..11].trim().parse()?;
        let x = coord(&line[30..38])?;
        let y = coord(&line[38..46])?;
        let z = coord(&line[46..54])?;

        atoms.push(Atom {
            chain,
            resseq,
            icode,
            resname,
            atom_name,
            serial,
            x,
            y,
            z,
        });
    }

    Ok(atoms)
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
            "H, HOH dropped; first altloc CB kept"
        );
        assert_eq!(atoms[2].chain, 'B');
        assert_eq!(atoms[2].x, 1.0, "altloc A coordinates win");
    }

    #[test]
    fn coordinates_round_through_f32() {
        let atoms = parse(SAMPLE).unwrap();
        // 16.967 has no exact f64 form; biopython sees its f32 rounding
        assert_eq!(atoms[1].x, f64::from(16.967_f32));
    }
}
