# ---------------------------------------------------------------------------
# Render HOMO and LUMO isosurfaces of lapachol from ORCA cube files.
#
# First generate the cubes with orca_plot (interactive):
#     orca_plot 02_sp_solvent.gbw -i
#       -> 1  (plot type: molecular orbital)
#       -> 2  (orbital number: use the HOMO index printed by orca_descriptors.py)
#       -> 5  (output format: Gaussian cube)
#       -> 4  (grid: 80 80 80)
#       -> 10 (generate)
#
# Then:  vmd -e render_orbitals.tcl
# ---------------------------------------------------------------------------

mol new lapachol_homo.cube type cube waitfor all
mol addfile lapachol_lumo.cube type cube waitfor all

# --- representation 0: the molecule itself, licorice
mol delrep 0 top
mol representation Licorice 0.12 24 24
mol color Name
mol selection {all}
mol material AOChalky
mol addrep top

# --- representations 1 and 2: HOMO, positive and negative lobes
foreach {iso colourid} {0.03 0 -0.03 1} {
    mol representation Isosurface $iso 0 0 0 1 1
    mol color ColorID $colourid
    mol selection {all}
    mol material Transparent
    mol addrep top
}

# --- representations 3 and 4: LUMO
foreach {iso colourid} {0.03 7 -0.03 3} {
    mol representation Isosurface $iso 1 0 0 1 1
    mol color ColorID $colourid
    mol selection {all}
    mol material Transparent
    mol addrep top
}

# --- presentation settings
display projection Orthographic
display depthcue off
display rendermode GLSL
axes location Off
color Display Background white

display resetview
rotate x by -60

# --- publication-quality output
render TachyonLOptiXInternal ../figures/orbitals.png

puts "Rendered to figures/orbitals.png"
puts "Toggle representations 1-2 (HOMO) against 3-4 (LUMO) in the GUI to compare."
