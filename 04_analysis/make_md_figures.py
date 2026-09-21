#!/usr/bin/env python3
"""
Generate the molecular-dynamics figures for the lectures and the README.

Produces four stand-alone panels, each sized for a slide:

    md_distance.png   lapachol–FAD minimum distance against time — the headline result
    md_contacts.png   contact occupancy of the site residues
    md_rg.png         radius of gyration per monomer — the diagnostic that proved
                      the protein was intact while the global RMSD said otherwise
    md_rmsd.png       backbone RMSD of chain A alone

Usage:
    python make_md_figures.py                         # light background
    python make_md_figures.py --dark                  # for dark slides
    python make_md_figures.py --top ../03_md/md_fix.pdb --traj ../03_md/md_fix.xtc
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import distances, rms

# Palette matching the lecture decks
INK = "#16212B"
CORAL = "#B5482A"
TEAL = "#1F6F63"
SAND = "#C98C6E"
MUTED = "#6B7280"

CONTACT_CUTOFF = 4.0   # angstrom, heavy atom to heavy atom
STACK_CUTOFF = 5.0     # angstrom, the line drawn on the distance plot


def style(dark: bool) -> dict:
    """Colours for a light or a dark slide background."""
    if dark:
        return dict(bg=INK, fg="#E8E4DC", grid="#2C3A47", muted="#8A9199")
    return dict(bg="#FFFFFF", fg=INK, grid="#E6E2DB", muted=MUTED)


def frame(ax, c):
    """Quiet the axes: no top or right spine, soft grid, muted ticks."""
    ax.set_facecolor(c["bg"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(c["grid"])
    ax.tick_params(colors=c["muted"], labelsize=11)
    ax.yaxis.label.set_color(c["fg"])
    ax.xaxis.label.set_color(c["fg"])
    ax.grid(axis="y", color=c["grid"], linewidth=0.8)
    ax.set_axisbelow(True)


def save(fig, path, c):
    fig.patch.set_facecolor(c["bg"])
    fig.tight_layout()
    fig.savefig(path, dpi=200, facecolor=c["bg"])
    plt.close(fig)
    print(f"  wrote {path}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--top", type=Path, default=Path("../03_md/md_fix.pdb"))
    ap.add_argument("--traj", type=Path, default=Path("../03_md/md_fix.xtc"))
    ap.add_argument("--out", type=Path, default=Path("../figures"))
    ap.add_argument("--dark", action="store_true", help="colours for a dark slide")
    args = ap.parse_args()

    c = style(args.dark)
    suffix = "_dark" if args.dark else ""
    args.out.mkdir(parents=True, exist_ok=True)

    u = mda.Universe(str(args.top), str(args.traj))
    protein = u.select_atoms("protein")
    n_half = len(protein) // 2
    chainA = protein[:n_half]
    chainB = protein[n_half:]
    lig = u.select_atoms("resname LIG")
    fad = u.select_atoms("resname FAD")

    print(f"Loaded {len(u.atoms)} atoms, {len(u.trajectory)} frames, "
          f"{u.trajectory.totaltime/1000:.1f} ns")
    print(f"  chain A {len(chainA)} atoms · chain B {len(chainB)} atoms · "
          f"LIG {len(lig)} · FAD {len(fad)}")

    # ---- single pass over the trajectory --------------------------------
    times, dmin, rgA, rgB, ncontacts = [], [], [], [], []
    tally: Counter = Counter()

    lig_heavy = lig.select_atoms("not name H*")
    prot_heavy = protein.select_atoms("not name H*")

    for ts in u.trajectory:
        times.append(ts.time / 1000.0)
        # minimum-image convention: immune to periodic-boundary wrapping
        d = distances.distance_array(lig.positions, fad.positions, box=ts.dimensions)
        dmin.append(float(d.min()))
        rgA.append(chainA.radius_of_gyration() / 10.0)
        rgB.append(chainB.radius_of_gyration() / 10.0)

        dp = distances.distance_array(lig_heavy.positions, prot_heavy.positions,
                                      box=ts.dimensions)
        near_idx = np.unique(np.where(dp.min(axis=0) < CONTACT_CUTOFF)[0])
        ncontacts.append(len(near_idx))
        for res in prot_heavy[near_idx].residues:
            tally[(int(res.resid), res.resname)] += 1

    times = np.array(times)
    dmin = np.array(dmin)
    rgA, rgB = np.array(rgA), np.array(rgB)
    nfr = len(u.trajectory)

    # ---- 1. lapachol–FAD minimum distance -------------------------------
    fig, ax = plt.subplots(figsize=(9, 4.6))
    frame(ax, c)
    ax.plot(times, dmin, linewidth=0.9, color=CORAL)
    ax.axhline(STACK_CUTOFF, color=TEAL, linestyle="--", linewidth=1.4)
    ax.text(times[-1], STACK_CUTOFF + 0.25, "van der Waals contact below 5 Å",
            ha="right", fontsize=11, color=TEAL)
    pct = 100.0 * (dmin < STACK_CUTOFF).sum() / nfr
    ax.text(0.02, 0.94, f"median {np.median(dmin):.2f} Å   ·   below 5 Å in {pct:.0f}% of frames",
            transform=ax.transAxes, fontsize=12, color=c["fg"], va="top")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Minimum distance (Å)")
    ax.set_title("Lapachol–FAD minimum distance", color=c["fg"], fontsize=14, loc="left")
    ax.set_ylim(0, max(12, dmin.max() * 1.1))
    save(fig, args.out / f"md_distance{suffix}.png", c)

    # ---- 2. contact occupancy -------------------------------------------
    top = tally.most_common(10)[::-1]
    labels = [f"{name}{resid}" for (resid, name), _ in top]
    values = [100.0 * n / nfr for _, n in top]
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    frame(ax, c)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", color=c["grid"], linewidth=0.8)
    ax.barh(labels, values, color=TEAL, height=0.68)
    for y, v in enumerate(values):
        ax.text(v + 1.2, y, f"{v:.0f}%", va="center", fontsize=11, color=c["muted"])
    ax.set_xlim(0, 100)
    ax.set_xlabel("Occupancy (% of frames)")
    ax.set_title(f"Residues within {CONTACT_CUTOFF:.0f} Å of lapachol",
                 color=c["fg"], fontsize=14, loc="left")
    save(fig, args.out / f"md_contacts{suffix}.png", c)

    # ---- 3. radius of gyration, the diagnostic --------------------------
    fig, ax = plt.subplots(figsize=(9, 4.6))
    frame(ax, c)
    ax.plot(times, rgA, linewidth=1.4, color=INK if not args.dark else "#E8E4DC",
            label="chain A")
    ax.plot(times, rgB, linewidth=1.4, color=SAND, label="chain B")
    ax.set_ylim(1.5, 2.5)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Radius of gyration (nm)")
    ax.set_title("Each monomer stays folded — the measure that needs no alignment",
                 color=c["fg"], fontsize=13.5, loc="left")
    leg = ax.legend(frameon=False, fontsize=11)
    for t in leg.get_texts():
        t.set_color(c["fg"])
    ax.text(0.02, 0.10, f"chain A {rgA.mean():.2f} ± {rgA.std():.2f} nm   ·   "
                        f"chain B {rgB.mean():.2f} ± {rgB.std():.2f} nm",
            transform=ax.transAxes, fontsize=11.5, color=c["muted"])
    save(fig, args.out / f"md_rg{suffix}.png", c)

    # ---- 4. backbone RMSD of chain A alone ------------------------------
    ca = chainA.select_atoms("name CA")
    sel = "index " + " ".join(map(str, ca.indices))
    analysis = rms.RMSD(u, select=sel, ref_frame=0).run()
    r = analysis.results.rmsd
    t_ns, bb = r[:, 1] / 1000.0, r[:, 2]
    half = bb[len(bb) // 2:]

    fig, ax = plt.subplots(figsize=(9, 4.6))
    frame(ax, c)
    ax.plot(t_ns, bb, linewidth=1.0, color=TEAL)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Backbone RMSD (Å)")
    ax.set_title("Chain A backbone RMSD", color=c["fg"], fontsize=14, loc="left")
    ax.set_ylim(0, max(4.0, bb.max() * 1.15))
    ax.text(0.02, 0.94, f"second half: {half.mean():.2f} ± {half.std():.2f} Å",
            transform=ax.transAxes, fontsize=12, color=c["fg"], va="top")
    save(fig, args.out / f"md_rmsd{suffix}.png", c)

    # ---- numbers for the README -----------------------------------------
    print("\nFor the README:")
    print(f"  lapachol–FAD minimum distance, median : {np.median(dmin):.2f} Å")
    print(f"  frames below 5 Å                      : {(dmin<5).sum()} of {nfr} ({pct:.0f}%)")
    print(f"  radius of gyration, chain A           : {rgA.mean():.2f} ± {rgA.std():.2f} nm")
    print(f"  backbone RMSD chain A, second half    : {half.mean():.2f} ± {half.std():.2f} Å")
    print(f"  mean heavy-atom contacts              : {np.mean(ncontacts):.1f}")


if __name__ == "__main__":
    main()
