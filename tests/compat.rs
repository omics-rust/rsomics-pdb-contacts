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
    pdb: &'static str,
    golden: &'static str,
    args: &'static [&'static str],
}

const CASES: &[Case] = &[
    Case {
        pdb: "tests/golden/1crn.pdb",
        golden: "tests/golden/1crn.residue_8.0.tsv",
        args: &["--cutoff", "8.0", "--level", "residue"],
    },
    Case {
        pdb: "tests/golden/1crn.pdb",
        golden: "tests/golden/1crn.atom_4.0.tsv",
        args: &["--cutoff", "4.0", "--level", "atom"],
    },
    Case {
        pdb: "tests/golden/1crn.pdb",
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
        pdb: "tests/golden/1crn.pdb",
        golden: "tests/golden/1crn.atom_ca_boundary.tsv",
        args: &[
            "--cutoff",
            "6.703862694776753",
            "--level",
            "atom",
            "--ca-only",
        ],
    },
    // Disordered atom: residue 1's CA has altloc A (occ 0.40, far) and B (occ
    // 0.60, near). biopython selects the higher-occupancy B, so it contacts
    // residues 2 and 3; the old first-seen-A parser saw none of those pairs.
    Case {
        pdb: "tests/golden/altloc.pdb",
        golden: "tests/golden/altloc.atom_ca_3.0.tsv",
        args: &["--cutoff", "3.0", "--level", "atom", "--ca-only"],
    },
];

#[test]
fn golden_pairs_match_biopython() {
    for case in CASES {
        let pdb = manifest(case.pdb);
        let pdb = pdb.to_str().unwrap();
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

/// Non-positive cutoffs must fail loudly the way biopython's
/// NeighborSearch.search_all raises `ValueError: Radius must be positive.`
#[test]
fn nonpositive_cutoff_errors_loud() {
    let pdb = manifest("tests/golden/1crn.pdb");
    let pdb = pdb.to_str().unwrap();
    for bad in ["0", "0.0", "-1", "-2.5"] {
        // `=` form so clap does not read a leading-minus value as a flag
        let cutoff = format!("--cutoff={bad}");
        let out = Command::new(bin())
            .args([pdb, &cutoff, "--level", "atom"])
            .output()
            .expect("spawn binary");
        assert!(
            !out.status.success(),
            "cutoff {bad} must exit non-zero, not silently proceed"
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("radius must be positive"),
            "cutoff {bad} stderr should name the positive-radius rule, got: {stderr}"
        );
    }
}

/// A valid but tiny cutoff explodes a naive uniform grid; the fallback must
/// return biopython's result (0 pairs on 1crn) without OOM/abort.
#[test]
fn tiny_cutoff_returns_empty_fast() {
    let pdb = manifest("tests/golden/1crn.pdb");
    let pdb = pdb.to_str().unwrap();
    for cut in ["0.001", "0.01", "0.02"] {
        for level in ["atom", "residue"] {
            let ours = run_ours(&[pdb, "--cutoff", cut, "--level", level]);
            let header = ours.lines().next().unwrap_or("");
            assert_eq!(
                ours.lines().count(),
                1,
                "cutoff {cut} level {level} must yield header only (0 pairs), got:\n{ours}"
            );
            assert!(header.starts_with("chainA"), "header preserved");
        }
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
    for case in CASES {
        let pdb = manifest(case.pdb);
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
