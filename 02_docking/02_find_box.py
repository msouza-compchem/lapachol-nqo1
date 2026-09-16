#!/usr/bin/env python3
"""
Determine the AutoDock Vina search box from the duroquinone position.

The receptor comes from 1D4A (apo human NQO1), but the substrate position is
only known from 1DXO (the duroquinone complex). Those two structures do not
share a coordinate frame - 1DXO contains two dimers in the asymmetric unit,
1D4A one - so 1DXO must be superposed onto 1D4A before its ligand coordinates
mean anything in the receptor's frame.

This script superposes the two on their common C-alpha atoms, applies the
resulting transformation to the whole mobile structure, and reports the
duroquinone centroid in the receptor frame, plus a box size that covers it
with margin.

Usage:
    python 02_find_box.py
"""

from __future__ import annotations

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align

REF_PDB = "1D4A.pdb"   # receptor frame: apo human NQO1
MOB_PDB = "1DXO.pdb"   # duroquinone complex, to be superposed
LIGAND = "DQN"         # duroquinone residue name
MARGIN = 8.0           # angstrom added around the ligand extent


def report_header(path: str) -> None:
    """Print the PDB TITLE/COMPND lines so the species is not assumed."""
    print(f"\n--- {path} ---")
    with open(path) as fh:
        for line in fh:
            if line.startswith(("TITLE", "COMPND", "SOURCE")):
                print("   ", line.rstrip())
            if line.startswith("ATOM"):
                break


def main() -> None:
    report_header(REF_PDB)
    report_header(MOB_PDB)

    ref = mda.Universe(REF_PDB)
    mob = mda.Universe(MOB_PDB)

    # Chains present in each structure
    ref_chains = sorted(set(ref.select_atoms("protein").chainIDs))
    mob_chains = sorted(set(mob.select_atoms("protein").chainIDs))
    print(f"\nChains in {REF_PDB}: {ref_chains}")
    print(f"Chains in {MOB_PDB}: {mob_chains}")

    # Superpose on the C-alpha atoms shared by chain A of both structures.
    # Restricting to one chain avoids any ambiguity about which copy of the
    # dimer corresponds to which.
    ref_ca = ref.select_atoms("protein and name CA and chainID A")
    mob_ca = mob.select_atoms("protein and name CA and chainID A")
    common = np.intersect1d(ref_ca.resids, mob_ca.resids)
    print(f"\nResidues common to chain A of both: {len(common)}")

    resid_list = " ".join(str(r) for r in common)
    sel = f"protein and name CA and chainID A and resid {resid_list}"

    old_rmsd, new_rmsd = align.alignto(mob, ref, select=sel)
    print(f"C-alpha RMSD before superposition: {old_rmsd:.2f} A")
    print(f"C-alpha RMSD after  superposition: {new_rmsd:.2f} A")
    if new_rmsd > 2.0:
        print("  WARNING: poor superposition. The two structures may not be the")
        print("  same protein - check the headers above before trusting the box.")

    # Duroquinone copies, now expressed in the receptor's coordinate frame
    print(f"\n{LIGAND} copies (receptor frame):")
    best = None
    for chain in sorted(set(mob.select_atoms(f"resname {LIGAND}").chainIDs)):
        lig = mob.select_atoms(f"resname {LIGAND} and chainID {chain}")
        c = lig.center_of_geometry()
        print(f"   chain {chain}: {len(lig)} atoms, "
              f"centre = {c[0]:8.3f} {c[1]:8.3f} {c[2]:8.3f}")
        if chain == "A":
            best = lig

    if best is None:
        raise SystemExit(f"No {LIGAND} found in chain A - inspect the structure.")

    centre = best.center_of_geometry()
    extent = best.positions.max(axis=0) - best.positions.min(axis=0)
    size = np.maximum(extent + 2 * MARGIN, 20.0)   # never smaller than 20 A

    print("\n" + "=" * 58)
    print("Paste into vina_config.txt:")
    print("=" * 58)
    print(f"center_x = {centre[0]:.3f}")
    print(f"center_y = {centre[1]:.3f}")
    print(f"center_z = {centre[2]:.3f}")
    print()
    print(f"size_x = {size[0]:.0f}")
    print(f"size_y = {size[1]:.0f}")
    print(f"size_z = {size[2]:.0f}")
    print("=" * 58)

    # Save the superposed structure so the pose can be inspected in VMD and
    # used later as the redocking reference.
    mob.atoms.write("1DXO_on_1D4A.pdb")
    best.write("duroquinone_ref.pdb")
    print("\nWrote 1DXO_on_1D4A.pdb and duroquinone_ref.pdb")
    print("duroquinone_ref.pdb is the crystal pose to compare your redocking against.")


if __name__ == "__main__":
    main()
