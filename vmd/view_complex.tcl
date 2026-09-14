# ---------------------------------------------------------------------------
# Visualise the NQO1-lapachol complex from the MD trajectory.
#     vmd -e view_complex.tcl
# ---------------------------------------------------------------------------

mol new ../03_md/md.tpr type tpr waitfor all
mol addfile ../03_md/md.xtc type xtc first 0 last -1 step 10 waitfor all

mol delrep 0 top

# --- protein: cartoon, coloured by secondary structure
mol representation NewCartoon 0.30 12 4.1 0
mol color Structure
mol selection {protein}
mol material AOChalky
mol addrep top

# --- FAD cofactor
mol representation Licorice 0.15 20 20
mol color ColorID 4
mol selection {resname FAD}
mol material AOChalky
mol addrep top

# --- lapachol
mol representation Licorice 0.20 24 24
mol color ColorID 1
mol selection {resname LIG}
mol material AOChalky
mol addrep top

# --- residues lining the site, within 4 A of the ligand
mol representation Licorice 0.10 16 16
mol color Name
mol selection {protein and same residue as (within 4 of resname LIG)}
mol material Opaque
mol addrep top

# --- the catalytic gate, Tyr128 and Phe232
mol representation Licorice 0.14 20 20
mol color ColorID 7
mol selection {protein and resid 128 232}
mol material Opaque
mol addrep top

display projection Orthographic
display depthcue off
axes location Off
color Display Background white

# Remove periodic-boundary jumps before judging anything visually
package require pbctools
pbc unwrap -all
pbc wrap -center com -centersel "protein" -compound residue -all

display resetview
puts "Loaded. Use the Graphics > Representations panel to isolate contacts."
