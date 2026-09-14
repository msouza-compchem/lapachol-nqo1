#!/usr/bin/env bash
# Build the solvated NQO1-lapachol complex for GROMACS.
#
# THE HARD PART, stated honestly: this system has TWO non-standard residues,
# the lapachol ligand and the FAD cofactor. Neither is in the stock CHARMM36
# force field. Plan for this before you start:
#   - lapachol -> CGenFF web server (cgenff.com), free for academic use.
#               Check the reported penalty scores. Anything above ~50 means the
#               analogy is poor and the parameters need manual attention.
#   - FAD      -> published CHARMM-compatible parameters exist; fetch them
#               rather than letting CGenFF guess at a cofactor this large.
# If parameterising FAD becomes a time sink, a defensible fallback is to run
# the MD on the apo-protein/ligand complex and say so explicitly in the README.
set -euo pipefail

# ---------------------------------------------------------------------------
# 1. Protein topology (download CHARMM36m into charmm36.ff/ first)
# ---------------------------------------------------------------------------
gmx pdb2gmx -f ../02_docking/receptor_raw.pdb -o protein.gro -water tip3p -ignh

# ---------------------------------------------------------------------------
# 2. Ligand topology
#    Upload best_pose.mol2 to cgenff.com, download the .str, then:
#      python cgenff_charmm2gmx.py LIG lig.mol2 lig.str charmm36.ff
# ---------------------------------------------------------------------------
# (produces lig.itp, lig.prm, lig_ini.pdb)

# ---------------------------------------------------------------------------
# 3. Merge protein + ligand into complex.gro and edit topol.top by hand:
#      #include "lig.prm"   right after the forcefield include
#      #include "lig.itp"   after the protein .itp includes
#      LIG   1              at the end of [ molecules ]
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 4. Box, solvate, neutralise at physiological ionic strength
# ---------------------------------------------------------------------------
gmx editconf -f complex.gro -o box.gro -c -d 1.0 -bt cubic
gmx solvate  -cp box.gro -cs spc216.gro -o solv.gro -p topol.top
gmx grompp   -f mdp/em.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top \
                      -pname NA -nname CL -neutral -conc 0.15

# ---------------------------------------------------------------------------
# 5. Minimisation and equilibration
# ---------------------------------------------------------------------------
gmx grompp -f mdp/em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun  -v -deffnm em

# Index group merging protein and ligand for temperature coupling
gmx make_ndx -f em.gro -o index.ndx << 'NDX'
1 | 13
q
NDX

gmx grompp -f mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun  -deffnm nvt

gmx grompp -f mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr
gmx mdrun  -deffnm npt

# ---------------------------------------------------------------------------
# 6. Production. Do NOT run this on the laptop.
# ---------------------------------------------------------------------------
gmx grompp -f mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr
echo
echo "Upload md.tpr to Colab and run there:"
echo "  gmx mdrun -deffnm md -nb gpu -bonded gpu -pme gpu -update gpu"
