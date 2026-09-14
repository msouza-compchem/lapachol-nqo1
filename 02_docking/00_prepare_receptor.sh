#!/usr/bin/env bash
# Prepare human NQO1 (PDB 1D4A) as a rigid receptor for AutoDock Vina.
#
# Design decisions worth defending in the README:
#   - FAD is KEPT. The isoalloxazine ring forms one wall of the substrate site;
#     removing it leaves a fictitious cavity and the docking becomes meaningless.
#   - Crystallographic waters are removed, except none are retained here.
#   - The biological unit is the HOMODIMER: the catalytic site sits at the
#     interface, so chains A and B are both required.
set -euo pipefail

# 1. Download the structures
wget -q https://files.rcsb.org/download/1D4A.pdb    # apo human NQO1, 1.7 A
wget -q https://files.rcsb.org/download/1DXO.pdb    # duroquinone complex (control)

# 2. Strip waters and other solvent, keep protein + FAD
grep -E "^(ATOM|HETATM)" 1D4A.pdb \
  | grep -v "HOH" \
  | awk '$5=="A" || $5=="B" || substr($0,22,1)=="A" || substr($0,22,1)=="B"' \
  > receptor_raw.pdb

# 3. Add hydrogens at pH 7.4 and convert to PDBQT
obabel receptor_raw.pdb -O receptor.pdbqt -xr -p 7.4

echo "Receptor written to receptor.pdbqt"
echo
echo "NEXT: determine the search box centre from the duroquinone position."
echo "Align 1DXO onto 1D4A in PyMOL or VMD, then read the ligand centroid:"
echo "  PyMOL:  load 1D4A.pdb; load 1DXO.pdb; align 1DXO, 1D4A"
echo "          print(cmd.centerofmass('1DXO and resn DQN'))"
echo "Put those coordinates into vina_config.txt."
