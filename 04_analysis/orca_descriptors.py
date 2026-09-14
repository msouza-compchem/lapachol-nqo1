#!/usr/bin/env python3
"""
Extract conceptual-DFT descriptors from ORCA output files.

Parses the frontier orbital energies and total energies produced by Stage 1 and
derives the quantities relevant to quinone reduction by NQO1:

    chemical potential      mu    = (E_HOMO + E_LUMO) / 2
    chemical hardness       eta   = (E_LUMO - E_HOMO) / 2
    electrophilicity index  omega = mu^2 / (2 * eta)
    adiabatic electron affinity   = E(neutral) - E(anion), both optimised

The electrophilicity index is the descriptor that tends to track the ability of
a quinone to accept a hydride, which is the step NQO1 catalyses.

Usage:
    python orca_descriptors.py [--dft-dir ../01_dft] [--out ../figures]
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

HARTREE_TO_EV = 27.211386245988


def parse_final_energy(path: Path) -> float:
    """Return the last 'FINAL SINGLE POINT ENERGY' in an ORCA output, in Hartree."""
    text = path.read_text(errors="ignore")
    matches = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
    if not matches:
        raise ValueError(f"No final energy found in {path}")
    return float(matches[-1])


def parse_frontier_orbitals(path: Path) -> tuple[float, float, int]:
    """
    Return (E_HOMO, E_LUMO, homo_index) in eV from the last ORBITAL ENERGIES
    block of an ORCA output.

    ORCA prints, per orbital:   NO   OCC   E(Eh)   E(eV)
    The HOMO is the last orbital with a non-zero occupation.
    """
    text = path.read_text(errors="ignore")
    blocks = text.split("ORBITAL ENERGIES")
    if len(blocks) < 2:
        raise ValueError(f"No orbital energies found in {path}")

    row = re.compile(r"^\s*(\d+)\s+([\d.]+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")
    orbitals: list[tuple[int, float, float]] = []
    for line in blocks[-1].splitlines():
        m = row.match(line)
        if m:
            idx, occ, _e_hartree, e_ev = m.groups()
            orbitals.append((int(idx), float(occ), float(e_ev)))

    if not orbitals:
        raise ValueError(f"Could not parse the orbital table in {path}")

    occupied = [o for o in orbitals if o[1] > 0.5]
    virtual = [o for o in orbitals if o[1] <= 0.5]
    if not occupied or not virtual:
        raise ValueError(f"Incomplete orbital table in {path}")

    homo = occupied[-1]
    lumo = virtual[0]
    return homo[2], lumo[2], homo[0]


def check_frequencies(path: Path) -> bool:
    """True if the optimisation reached a genuine minimum (no imaginary modes)."""
    if not path.exists():
        return False
    text = path.read_text(errors="ignore")
    if "VIBRATIONAL FREQUENCIES" not in text:
        return False
    return "***imaginary mode***" not in text


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dft-dir", default="../01_dft", type=Path)
    ap.add_argument("--out", default="../figures", type=Path)
    args = ap.parse_args()

    dft = args.dft_dir
    args.out.mkdir(parents=True, exist_ok=True)

    opt_out = dft / "01_opt_freq.out"
    sp_out = dft / "02_sp_solvent.out"
    anion_out = dft / "03_anion_opt.out"

    # --- sanity check first: an imaginary frequency invalidates everything after
    is_minimum = check_frequencies(opt_out)
    print(f"Optimised structure is a true minimum: {is_minimum}")
    if not is_minimum:
        print("  -> Do not report these descriptors until the geometry is fixed.")

    # --- frontier orbitals from the solvated single point
    e_homo, e_lumo, homo_idx = parse_frontier_orbitals(sp_out)
    gap = e_lumo - e_homo
    mu = 0.5 * (e_homo + e_lumo)
    eta = 0.5 * (e_lumo - e_homo)
    omega = mu**2 / (2 * eta)

    print(f"\nHOMO is orbital #{homo_idx} (use {homo_idx}/{homo_idx + 1} in the "
          f"%plots block to render the cubes)")

    # --- adiabatic electron affinity
    ea = None
    if anion_out.exists():
        e_neutral = parse_final_energy(sp_out)
        e_anion = parse_final_energy(anion_out)
        ea = (e_neutral - e_anion) * HARTREE_TO_EV

    results = {
        "E_HOMO (eV)": e_homo,
        "E_LUMO (eV)": e_lumo,
        "HOMO-LUMO gap (eV)": gap,
        "Chemical potential mu (eV)": mu,
        "Chemical hardness eta (eV)": eta,
        "Electrophilicity index omega (eV)": omega,
    }
    if ea is not None:
        results["Adiabatic electron affinity (eV)"] = ea

    df = pd.DataFrame.from_dict(results, orient="index", columns=["Value"])
    df["Value"] = df["Value"].round(3)
    print("\n" + df.to_string())

    df.to_csv(args.out / "dft_descriptors.csv")

    # --- orbital energy diagram
    fig, ax = plt.subplots(figsize=(3.2, 4.6))
    for energy, label, colour in [
        (e_homo, f"HOMO\n{e_homo:.2f} eV", "#1f77b4"),
        (e_lumo, f"LUMO\n{e_lumo:.2f} eV", "#d62728"),
    ]:
        ax.hlines(energy, 0.25, 0.75, colors=colour, linewidth=3)
        ax.text(0.80, energy, label, va="center", fontsize=9, color=colour)

    ax.annotate(
        "", xy=(0.5, e_lumo), xytext=(0.5, e_homo),
        arrowprops=dict(arrowstyle="<->", color="grey", linewidth=1.2),
    )
    ax.text(0.44, mu, f"gap\n{gap:.2f} eV", ha="right", va="center",
            fontsize=9, color="grey")

    ax.set_xlim(0, 1.6)
    ax.set_ylim(e_homo - 1.5, e_lumo + 1.5)
    ax.set_ylabel("Energy (eV)")
    ax.set_xticks([])
    ax.set_title("Lapachol frontier orbitals\nB3LYP-D4/def2-TZVP, CPCM(water)",
                 fontsize=10)
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.out / "frontier_orbitals.png", dpi=300)
    print(f"\nWrote {args.out / 'dft_descriptors.csv'} and "
          f"{args.out / 'frontier_orbitals.png'}")


if __name__ == "__main__":
    main()
