use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_pdb_contacts::contacts::{Level, contacts};
use rsomics_pdb_contacts::pdb;

fn load() -> Vec<pdb::Atom> {
    let path =
        std::env::var("CONTACTS_BENCH_PDB").unwrap_or_else(|_| "tests/golden/1crn.pdb".to_string());
    let text = std::fs::read_to_string(path).expect("read bench pdb");
    pdb::parse(&text).expect("parse")
}

fn bench(c: &mut Criterion) {
    let atoms = load();
    c.bench_function("residue_8.0", |b| {
        b.iter(|| contacts(black_box(&atoms), 8.0, Level::Residue, false, 0).unwrap());
    });
    c.bench_function("atom_8.0", |b| {
        b.iter(|| contacts(black_box(&atoms), 8.0, Level::Atom, false, 0).unwrap());
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
