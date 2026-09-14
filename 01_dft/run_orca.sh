#!/usr/bin/env bash
# Run the full DFT stage. Adjust ORCA_BIN to your installation.
set -euo pipefail

ORCA_BIN="${ORCA_BIN:-$(command -v orca)}"
cd "$(dirname "$0")"

echo "==> [0/3] fetching structure"
bash 00_get_structure.sh

echo "==> [1/3] optimisation + frequencies"
"$ORCA_BIN" 01_opt_freq.inp > 01_opt_freq.out

# ORCA writes the optimised geometry to <basename>.xyz
cp 01_opt_freq.xyz ligand_opt.xyz

echo "    checking for imaginary frequencies..."
grep -A5 "VIBRATIONAL FREQUENCIES" 01_opt_freq.out | head -20
if grep -q "\*\*\*imaginary mode\*\*\*" 01_opt_freq.out; then
  echo "    WARNING: imaginary frequency found - this is NOT a minimum."
  echo "    Displace along the imaginary mode and re-optimise before continuing."
  exit 1
fi

echo "==> [2/3] single point, def2-TZVP + CPCM(water)"
"$ORCA_BIN" 02_sp_solvent.inp > 02_sp_solvent.out

echo "==> [3/3] radical anion optimisation"
"$ORCA_BIN" 03_anion_opt.inp > 03_anion_opt.out

echo "Done. Now run: python ../04_analysis/orca_descriptors.py"
