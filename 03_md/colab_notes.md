# Running the production MD on a free cloud GPU

The laptop handles everything except the 100 ns production run. On 4 CPU cores a
~50,000-atom system delivers roughly 1-3 ns/day; a single free-tier T4 delivers
30-80 ns/day.

## Workflow

1. Build, minimise and equilibrate locally (`build_system.sh`). This is all cheap.
2. Produce `md.tpr` locally. A `.tpr` is self-contained: it carries the topology,
   the parameters and the coordinates, so it is the only file you need to upload.
3. Upload `md.tpr` to Google Drive, mount the Drive in Colab, and run:

```python
!gmx mdrun -deffnm md -nb gpu -bonded gpu -pme gpu -update gpu -nstlist 100
```

4. Colab sessions time out. Checkpoint and restart rather than losing work:

```python
!gmx mdrun -deffnm md -cpi md.cpt -append -maxh 11 \
           -nb gpu -bonded gpu -pme gpu -update gpu
```

Write checkpoints straight to the mounted Drive so a disconnect costs minutes,
not days.

5. Download `md.xtc` and `md.tpr` and analyse locally. A 100 ns trajectory saved
   every 10 ps with only protein and ligand written out is a few hundred MB —
   comfortable on the laptop, and not committed to git.

## Ready-made notebooks

The **"Making it Rain"** notebook collection (Arantes, Polêto, Pedebos &
Ligabue-Braun) provides cloud-based GROMACS and OpenMM simulation notebooks built
for exactly this situation. Worth reading even if you write your own cells.
