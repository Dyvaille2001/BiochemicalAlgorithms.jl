# 1. Load the Minimized PDBs as reference structures
#and immediately strip hydrogens to match the "heavy" selection

load 5TGZ_minimized_cb1.pdb, wt_ref
remove (wt_ref and hydro)

load 5TGZ_F238L_repacked_minimized_cb1.pdb, mut_ref
remove (mut_ref and hydro)

# 2. Load the Mode 7 trajectories onto the structures
# This assumes the .xyz files have the same number of atoms as the PDBs
load_traj 5TGZ_heavy_mode_7.xyz, wt_ref
load_traj 5TGZ_F238L_repacked_heavy_mode_7.xyz, mut_ref

# 3. Clean up the view
hide everything
show cartoon, wt_ref or mut_ref
color gray80, wt_ref
color slate, mut_ref

# 4. Create the Porcupine Plot (Vectors)
# We use the 'modevectors' tool if installed, or simple arrows.
# Below creates a simple visualization of the movement
set cartoon_trace_atoms, 1
preset.technical(selection='all')

# 5. Side-by-side view
grid_mode 1
orient

# 6. Highlight the mutation site
show sticks, (mut_ref and resi 237)
color yellow, (mut_ref and resi 237)