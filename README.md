# rsomics-pdb-contacts

Residue/atom contact pairs within a distance cutoff from a PDB, value-exact to
`Bio.PDB.NeighborSearch.search_all`.

```bash
rsomics-pdb-contacts structure.pdb --cutoff 8.0 --level residue
rsomics-pdb-contacts structure.pdb --cutoff 5.0 --level atom
rsomics-pdb-contacts structure.pdb --cutoff 8.0 --ca-only --min-seq-sep 4
rsomics-pdb-contacts structure.pdb --cutoff 8.0 --json
```

## What a contact is

A pair of entities is a contact when at least one atom pair is within `--cutoff`
Angstrom. The boundary is **inclusive**: a pair whose distance is exactly equal
to the cutoff is included (biopython compares squared distance `r <= radius**2`).
Self-pairs never occur, and at residue level the atom pairs inside a single
residue are dropped.

- `--level residue` (default): residue pairs that have any atom pair within the
  cutoff, deduplicated.
- `--level atom`: the atom pairs themselves.
- `--ca-only`: restrict the atom universe to Cα atoms (the Cα-Cα convention).
- `--min-seq-sep K`: skip intra-chain pairs whose residue-sequence indices differ
  by less than `K` (0 = no filter).
- `--metric euclidean`: the only defined metric.

The atom universe matches biopython's default `PDBParser` feeding
`NeighborSearch`: first model, first altloc per residue, standard amino-acid
residues, hydrogens dropped. Coordinates are parsed through 32-bit floats, as
biopython stores them, so the inclusive boundary is bit-exact.

Output is a TSV with a fixed header, sorted by `(chainA, resA, chainB, resB)`
(atom level additionally by atom serial). `--json` emits the single
rsomics-common envelope.

## Origin

This crate is an independent Rust reimplementation of the contact search in
`Bio.PDB.NeighborSearch.search_all` based on:

- The biopython `NeighborSearch` / `kdtrees` behaviour (BSD-3-Clause), in
  particular the inclusive `r <= radius**2` neighbor predicate, the
  `index1 < index2` atom-pair canonicalisation, and the `level='R'`
  parent-pair folding (`_get_unique_parent_pairs`).
- Black-box behaviour testing against biopython 1.87.

The spatial index is a uniform cell-list, not biopython's KD-tree; only the
contact set and its ordering are required to match.

License: MIT OR Apache-2.0.
Upstream credit: Biopython (https://biopython.org/, BSD-3-Clause).
