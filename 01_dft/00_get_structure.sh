#!/usr/bin/env bash
# Fetch the 3D structure of lapachol from PubChem and convert to XYZ.
# PubChem's 3D record is a reasonable starting geometry; DFT will refine it.
set -euo pipefail

NAME="lapachol"

curl -L -o ligand_pubchem.sdf \
  "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${NAME}/SDF?record_type=3d"

# Add hydrogens explicitly and generate a clean XYZ for ORCA
obabel ligand_pubchem.sdf -O ligand.xyz -h --gen3d

echo "Wrote ligand.xyz"
echo "Heavy-atom count:"
grep -vc "^ *H" ligand.xyz || true

# OPTIONAL but recommended: conformational search with CREST (GFN2-xTB).
# Lapachol has a flexible prenyl chain; starting DFT from the wrong conformer
# is the single most common way to get a meaningless HOMO/LUMO.
#   crest ligand.xyz --gfn2 --alpb water -T 4
# Then use crest_best.xyz as the DFT starting point.
