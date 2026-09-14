#!/usr/bin/env bash
# Full docking stage, including protocol validation by redocking.
set -euo pipefail
cd "$(dirname "$0")"

bash 00_prepare_receptor.sh
bash 01_prepare_ligand.sh

echo "==> [1/2] control: redocking duroquinone"
# A protocol that cannot reproduce a known crystal pose cannot be trusted
# for an unknown one. Target: RMSD below 2.0 A against the 1DXO pose.
vina --config vina_config.txt \
     --ligand duroquinone.pdbqt \
     --out    duroquinone_redock.pdbqt \
     > duroquinone_redock.log

echo "==> [2/2] docking lapachol"
vina --config vina_config.txt \
     --out lapachol_poses.pdbqt \
     > lapachol_docking.log

echo
echo "Best scores:"
grep -A12 "mode |" lapachol_docking.log | head -15
echo
echo "Extract the best pose for MD:"
echo "  obabel lapachol_poses.pdbqt -O best_pose.pdb -f 1 -l 1"
