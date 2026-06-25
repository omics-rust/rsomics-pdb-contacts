//! Differential compat against Bio.PDB.NeighborSearch.search_all.
//!
//! The committed goldens (tests/golden/1crn.*.tsv) were captured from biopython
//! 1.87 via tests/contacts_oracle.py. The always-run tests reproduce them from
//! the rsomics binary with no biopython on the machine — they are the
//! authoritative gate. The live-oracle test re-derives the golden when biopython
//! is on PATH (loud-skip otherwise) and asserts byte equality, so a future
//! biopython that changed the contact definition would be caught.
//!
//! Contact set is value-exact: the SET of pairs and the row order match
//! biopython exactly. The cutoff boundary is inclusive (distance <= cutoff);
//! the boundary golden pins a pair whose distance equals the cutoff exactly.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-pdb-contacts"))
}

fn manifest(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn run_ours(args: &[&str]) -> String {
    let out = Command::new(bin())
        .args(args)
        .output()
        .expect("spawn binary");
    assert!(
        out.status.success(),
        "binary failed ({:?}): {}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("utf8 stdout")
}

struct Case {
    golden: &'static str,
    args: &'static [&'static str],
}

const CASES: &[Case] = &[
    Case {
        golden: "tests/golden/1crn.residue_8.0.tsv",
        args: &["--cutoff", "8.0", "--level", "residue"],
    },
    Case {
        golden: "tests/golden/1crn.atom_4.0.tsv",
        args: &["--cutoff", "4.0", "--level", "atom"],
    },
    Case {
        golden: "tests/golden/1crn.residue_ca_8.0_sep2.tsv",
        args: &[
            "--cutoff",
            "8.0",
            "--level",
            "residue",
            "--ca-only",
            "--min-seq-sep=2",
        ],
    },
    Case {
        golden: "tests/golden/1crn.atom_ca_boundary.tsv",
        args: &[
            "--cutoff",
            "6.703862694776753",
            "--level",
            "atom",
            "--ca-only",
        ],
    },
];

#[test]
fn golden_pairs_match_biopython() {
    let pdb = manifest("tests/golden/1crn.pdb");
    let pdb = pdb.to_str().unwrap();
    for case in CASES {
        let mut args: Vec<&str> = vec![pdb];
        args.extend_from_slice(case.args);
        let ours = run_ours(&args);
        let golden = std::fs::read_to_string(manifest(case.golden)).expect("read golden");
        assert_eq!(
            ours, golden,
            "contact-pair set/order mismatch for {:?}",
            case.args
        );
    }
}

#[test]
fn boundary_is_inclusive() {
    let pdb = manifest("tests/golden/1crn.pdb");
    let pdb = pdb.to_str().unwrap();
    let cut = "6.703862694776753";
    let included = run_ours(&[pdb, "--cutoff", cut, "--level", "atom", "--ca-only"]);
    // one ULP below the exact pair distance must drop exactly that pair
    let excluded = run_ours(&[
        pdb,
        "--cutoff",
        "6.703862694776752",
        "--level",
        "atom",
        "--ca-only",
    ]);
    let n_in = included.lines().count();
    let n_out = excluded.lines().count();
    assert_eq!(
        n_in,
        n_out + 1,
        "exact-distance pair must be included at cutoff == distance"
    );
}

#[test]
fn live_oracle_biopython() {
    let Some(py) = which_biopython() else {
        eprintln!("SKIP live_oracle_biopython: no biopython on PATH");
        return;
    };
    let oracle = manifest("tests/contacts_oracle.py");
    let pdb = manifest("tests/golden/1crn.pdb");
    for case in CASES {
        let mut oargs: Vec<String> = vec![
            oracle.to_string_lossy().into_owned(),
            pdb.to_string_lossy().into_owned(),
        ];
        oargs.push(cutoff_of(case.args));
        oargs.push(level_of(case.args));
        for extra in extras_of(case.args) {
            oargs.push(extra);
        }
        let out = Command::new(&py)
            .args(&oargs)
            .output()
            .expect("spawn oracle");
        assert!(
            out.status.success(),
            "oracle failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let expected = String::from_utf8(out.stdout).expect("utf8");

        let mut args: Vec<&str> = vec![pdb.to_str().unwrap()];
        args.extend_from_slice(case.args);
        let ours = run_ours(&args);
        assert_eq!(ours, expected, "live oracle mismatch for {:?}", case.args);
    }
}

fn cutoff_of(args: &[&str]) -> String {
    let i = args.iter().position(|a| *a == "--cutoff").unwrap();
    args[i + 1].to_string()
}

fn level_of(args: &[&str]) -> String {
    let i = args.iter().position(|a| *a == "--level").unwrap();
    args[i + 1].to_string()
}

fn extras_of(args: &[&str]) -> Vec<String> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < args.len() {
        match args[i] {
            "--cutoff" | "--level" => i += 2,
            other => {
                out.push(other.to_string());
                i += 1;
            }
        }
    }
    out
}

fn which_biopython() -> Option<String> {
    for py in [
        "/opt/homebrew/Caskroom/miniforge/base/envs/rs-up/bin/python",
        "python3",
        "python",
    ] {
        let ok = Command::new(py)
            .args(["-c", "import Bio.PDB.NeighborSearch"])
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);
        if ok {
            return Some(py.to_string());
        }
    }
    None
}
