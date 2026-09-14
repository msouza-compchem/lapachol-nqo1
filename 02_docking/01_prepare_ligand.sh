#!/usr/bin/env bash
# Convert the DFT-optimised lapachol geometry into a Vina-ready PDBQT.
# Using the DFT geometry rather than a raw PubChem structure is a small but
# real methodological improvement, and it links Stage 1 to Stage 2.
set -euo pipefail

obabel ../01_dft/ligand_opt.xyz -O ligand.pdbqt -p 7.4 --partialcharge gasteiger

# Control ligand for protocol validation
obabel 1DXO.pdb -O duroquinone.pdbqt -p 7.4 --partialcharge gasteiger 2>/dev/null || \
  echo "Extract the DQN residue from 1DXO.pdb manually before redocking."

echo "Ligands prepared."
