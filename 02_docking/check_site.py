import MDAnalysis as mda
import numpy as np

u = mda.Universe('receptor_raw.pdb')

for name, c in [("cadeia A", [19.754, -6.423, 5.313]),
                ("cadeia C", [-1.339, 14.682, 5.149])]:
    near = u.select_atoms(f'point {c[0]} {c[1]} {c[2]} 10')
    fad = u.select_atoms(f'resname FAD and point {c[0]} {c[1]} {c[2]} 12')
    print(f"\n--- sitio da {name}: {c} ---")
    print(f"  atomos a 10 A   : {len(near)}")
    print(f"  atomos de FAD   : {len(fad)}")
    if len(near):
        res = sorted({(r.resname, r.resid) for r in near.residues},
                     key=lambda x: x[1])
        print(f"  residuos ({len(res)}): {res[:12]}")
