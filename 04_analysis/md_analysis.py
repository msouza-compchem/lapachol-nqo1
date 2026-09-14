#!/usr/bin/env python3
"""
Analyse the NQO1-lapachol molecular dynamics trajectory.

Produces the four panels that answer the only question the MD is there to
answer: does the docked pose survive, and what holds it in place?

    1. Backbone RMSD          - did the protein itself stay folded and converge?
    2. Ligand RMSD            - did the pose drift out of the site?
    3. Per-residue RMSF       - which parts of the protein are mobile
                                (watch Tyr128 and Phe232, the catalytic gate)
    4. Protein-ligand contacts - how many heavy-atom contacts persist, and which
                                residues are responsible

Usage:
    python md_analysis.py --topology ../03_md/md.tpr --trajectory ../03_md/md.xtc
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms

# Residues of mechanistic interest in human NQO1 (Faig et al., PNAS 2000)
GATE_RESIDUES = {128: "Tyr128", 232: "Phe232"}

CONTACT_CUTOFF = 4.0  # angstrom, heavy-atom to heavy-atom


def load_universe(topology: Path, trajectory: Path) -> mda.Universe:
    u = mda.Universe(str(topology), str(trajectory))
    print(f"Loaded {len(u.atoms)} atoms, {len(u.trajectory)} frames, "
          f"{u.trajectory.totaltime / 1000:.1f} ns")
    return u


def backbone_and_ligand_rmsd(u: mda.Universe, ligand_sel: str) -> pd.DataFrame:
    """RMSD of the protein backbone and of the ligand, after aligning on the protein."""
    analysis = rms.RMSD(
        u,
        select="backbone",
        groupselections=[ligand_sel],
        ref_frame=0,
    ).run()
    data = analysis.results.rmsd
    return pd.DataFrame({
        "time_ns": data[:, 1] / 1000.0,
        "backbone_rmsd_A": data[:, 2],
        "ligand_rmsd_A": data[:, 3],
    })


def per_residue_rmsf(u: mda.Universe) -> pd.DataFrame:
    """RMSF of C-alpha atoms, computed on the trajectory aligned to its own average."""
    protein = u.select_atoms("protein")
    aligner = align.AlignTraj(u, u, select="protein and name CA", in_memory=True)
    aligner.run()
    calphas = protein.select_atoms("name CA")
    rmsf = rms.RMSF(calphas).run()
    return pd.DataFrame({
        "resid": calphas.resids,
        "resname": calphas.resnames,
        "rmsf_A": rmsf.results.rmsf,
    })


def contact_analysis(u: mda.Universe, ligand_sel: str) -> tuple[pd.DataFrame, Counter]:
    """Count heavy-atom contacts per frame and tally which residues contribute."""
    ligand = u.select_atoms(f"({ligand_sel}) and not name H*")
    counts: list[int] = []
    times: list[float] = []
    residue_tally: Counter = Counter()

    for ts in u.trajectory:
        near = u.select_atoms(
            f"protein and not name H* and around {CONTACT_CUTOFF} group lig",
            lig=ligand,
            updating=False,
        )
        counts.append(len(near))
        times.append(ts.time / 1000.0)
        for res in near.residues:
            residue_tally[(int(res.resid), res.resname)] += 1

    n_frames = len(u.trajectory)
    # convert raw tallies into occupancy as a fraction of the trajectory
    occupancy = Counter({k: v / n_frames for k, v in residue_tally.items()})
    return pd.DataFrame({"time_ns": times, "n_contacts": counts}), occupancy


def plot_all(rmsd: pd.DataFrame, rmsf: pd.DataFrame,
             contacts: pd.DataFrame, occupancy: Counter, outdir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))

    # --- panel 1: backbone RMSD
    ax = axes[0, 0]
    ax.plot(rmsd["time_ns"], rmsd["backbone_rmsd_A"], linewidth=0.8, color="#1f77b4")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Backbone RMSD (Å)")
    ax.set_title("Protein stability")

    # --- panel 2: ligand RMSD
    ax = axes[0, 1]
    ax.plot(rmsd["time_ns"], rmsd["ligand_rmsd_A"], linewidth=0.8, color="#d62728")
    ax.axhline(2.0, color="grey", linestyle="--", linewidth=0.8)
    ax.text(rmsd["time_ns"].iloc[-1], 2.05, "pose retained below 2 Å",
            ha="right", fontsize=8, color="grey")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Ligand RMSD (Å)")
    ax.set_title("Docked pose stability")

    # --- panel 3: RMSF with the catalytic gate highlighted
    ax = axes[1, 0]
    ax.plot(rmsf["resid"], rmsf["rmsf_A"], linewidth=0.8, color="#2ca02c")
    for resid, label in GATE_RESIDUES.items():
        if resid in set(rmsf["resid"]):
            value = float(rmsf.loc[rmsf["resid"] == resid, "rmsf_A"].iloc[0])
            ax.plot(resid, value, "o", color="#d62728", markersize=5)
            ax.annotate(label, (resid, value), textcoords="offset points",
                        xytext=(4, 6), fontsize=8, color="#d62728")
    ax.set_xlabel("Residue number")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title("Per-residue flexibility")

    # --- panel 4: contact occupancy
    ax = axes[1, 1]
    top = occupancy.most_common(10)[::-1]
    labels = [f"{name}{resid}" for (resid, name), _ in top]
    values = [frac for _, frac in top]
    ax.barh(labels, values, color="#9467bd")
    ax.set_xlabel("Occupancy (fraction of trajectory)")
    ax.set_xlim(0, 1)
    ax.set_title(f"Contacts within {CONTACT_CUTOFF} Å")

    fig.suptitle("NQO1–lapachol, 100 ns, CHARMM36m/TIP3P, 300 K", fontsize=11)
    fig.tight_layout()
    fig.savefig(outdir / "md_summary.png", dpi=300)
    print(f"Wrote {outdir / 'md_summary.png'}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--topology", type=Path, default=Path("../03_md/md.tpr"))
    ap.add_argument("--trajectory", type=Path, default=Path("../03_md/md.xtc"))
    ap.add_argument("--ligand-sel", default="resname LIG",
                    help="MDAnalysis selection string for the ligand")
    ap.add_argument("--out", type=Path, default=Path("../figures"))
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    u = load_universe(args.topology, args.trajectory)

    print("Computing RMSD...")
    rmsd = backbone_and_ligand_rmsd(u, args.ligand_sel)

    print("Computing RMSF...")
    rmsf = per_residue_rmsf(u)

    print("Computing contacts...")
    contacts, occupancy = contact_analysis(u, args.ligand_sel)

    # --- numbers for the README table: average over the last half of the run
    second_half = rmsd[rmsd["time_ns"] >= rmsd["time_ns"].max() / 2]
    print("\nEquilibrated averages (second half of the trajectory):")
    print(f"  backbone RMSD : {second_half['backbone_rmsd_A'].mean():.2f} "
          f"+/- {second_half['backbone_rmsd_A'].std():.2f} A")
    print(f"  ligand RMSD   : {second_half['ligand_rmsd_A'].mean():.2f} "
          f"+/- {second_half['ligand_rmsd_A'].std():.2f} A")
    print(f"  mean contacts : {contacts['n_contacts'].mean():.1f}")

    print("\nMost persistent contacts:")
    for (resid, name), frac in occupancy.most_common(10):
        print(f"  {name}{resid:<5} {frac:5.1%}")

    rmsd.to_csv(args.out / "rmsd.csv", index=False)
    rmsf.to_csv(args.out / "rmsf.csv", index=False)
    contacts.to_csv(args.out / "contacts.csv", index=False)
    plot_all(rmsd, rmsf, contacts, occupancy, args.out)


if __name__ == "__main__":
    main()
