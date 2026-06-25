"""Reference contact-pair generator using Bio.PDB.NeighborSearch.

Emits the same sorted contact-pair table rsomics-pdb-contacts produces, so the
committed golden and the live-oracle test share one source of truth.

Contact definition (Bio.PDB.NeighborSearch.search_all, biopython 1.87):
  - all-atom search: every atom of every standard residue is a query point.
  - a pair (i, j) is a contact iff squared distance <= cutoff**2
    (kdtrees.c KDTree_test_neighbors: `r <= self->_neighbor_radius_sq`),
    so the cutoff boundary is INCLUSIVE and self-pairs (i == j) never occur
    (inner loop starts at j = i + 1).
  - level='A' returns atom pairs (index1 < index2 canonical order).
  - level='R' folds atom pairs to residue pairs via _get_unique_parent_pairs:
    same-residue pairs are dropped, the smaller residue (by (model, chain,
    (hetfield, resseq, icode))) comes first, and duplicates collapse to a set.

This script reproduces the same residue/atom universe rsomics uses: standard
amino-acid residues only, first model, first altloc, no hydrogens.
"""

import sys
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch

STD = set(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE "
    "LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
)


def is_hydrogen(atom):
    el = (atom.element or "").strip()
    if el:
        return el == "H" or el == "D"
    name = atom.get_name().strip()
    return name[:1] in ("H", "D")


def collect_atoms(pdb_path, ca_only):
    parser = PDBParser(QUIET=True)
    model = parser.get_structure("x", pdb_path)[0]
    atoms = []
    for chain in model:
        for res in chain:
            if res.get_resname() not in STD:
                continue
            for atom in res:
                if is_hydrogen(atom):
                    continue
                if ca_only and atom.get_name() != "CA":
                    continue
                atoms.append(atom)
    return atoms


def res_key(residue):
    chain = residue.get_parent().id
    het, seq, icode = residue.id
    return (chain, het, seq, icode)


def res_label(residue):
    chain = residue.get_parent().id
    het, seq, icode = residue.id
    ic = icode.strip()
    seqtag = f"{seq}{ic}" if ic else f"{seq}"
    return f"{chain}\t{seqtag}\t{residue.get_resname()}"


def atom_key(atom):
    res = atom.get_parent()
    chain, het, seq, icode = res_key(res)
    return (chain, het, seq, icode, atom.get_serial_number())


def atom_label(atom):
    res = atom.get_parent()
    chain = res.get_parent().id
    het, seq, icode = res.id
    ic = icode.strip()
    seqtag = f"{seq}{ic}" if ic else f"{seq}"
    return f"{chain}\t{seqtag}\t{res.get_resname()}\t{atom.get_name()}"


def seq_index(residue):
    return residue.id[1]


def main():
    pdb_path = sys.argv[1]
    cutoff = float(sys.argv[2])
    level = sys.argv[3] if len(sys.argv) > 3 else "residue"
    ca_only = "--ca-only" in sys.argv
    min_seq_sep = 0
    for tok in sys.argv:
        if tok.startswith("--min-seq-sep="):
            min_seq_sep = int(tok.split("=", 1)[1])

    atoms = collect_atoms(pdb_path, ca_only)
    ns = NeighborSearch(atoms)

    if level == "atom":
        pairs = ns.search_all(cutoff, level="A")
        rows = set()
        for a, b in pairs:
            ra, rb = a.get_parent(), b.get_parent()
            if min_seq_sep > 0 and ra.get_parent().id == rb.get_parent().id:
                if abs(seq_index(ra) - seq_index(rb)) < min_seq_sep:
                    continue
            ka, kb = atom_key(a), atom_key(b)
            if kb < ka:
                a, b, ka, kb = b, a, kb, ka
            rows.add((ka, kb, atom_label(a), atom_label(b)))
        out = ["\t".join([la, lb]) for (_, _, la, lb) in sorted(rows)]
        print("chainA\tresA\tresnameA\tatomA\tchainB\tresB\tresnameB\tatomB")
    else:
        pairs = ns.search_all(cutoff, level="R")
        rows = set()
        for ra, rb in pairs:
            if min_seq_sep > 0 and ra.get_parent().id == rb.get_parent().id:
                if abs(seq_index(ra) - seq_index(rb)) < min_seq_sep:
                    continue
            ka, kb = res_key(ra), res_key(rb)
            if kb < ka:
                ra, rb, ka, kb = rb, ra, kb, ka
            rows.add((ka, kb, res_label(ra), res_label(rb)))
        out = ["\t".join([la, lb]) for (_, _, la, lb) in sorted(rows)]
        print("chainA\tresA\tresnameA\tchainB\tresB\tresnameB")

    for line in out:
        print(line)


if __name__ == "__main__":
    main()
