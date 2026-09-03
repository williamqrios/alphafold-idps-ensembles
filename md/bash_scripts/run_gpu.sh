#!/bin/bash
OUTDIR="out"
cd $OUTDIR
gmx_gpu grompp -f ../em.mdp -c ../protein_neutral.gro -p ../topol.top -o em.tpr
gmx_gpu mdrun -v -deffnm em

# Equilibration NVT for 1 ns (full restraint)
gmx_gpu grompp -f ../nvt.mdp -c em.gro -r em.gro -p ../topol.top -o nvt.tpr
gmx_gpu mdrun -v -deffnm nvt

# Equilibration NPT for 1 ns (full restraint), -maxwarn required due to Berendsen barostat
if [ -f nvt.xtc ]; then
    gmx_gpu grompp -f ../npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p ../topol.top -o npt.tpr -maxwarn 2
    gmx_gpu mdrun -v -deffnm npt
fi

# Full run 
if [ -f npt.xtc ]; then
   gmx_gpu grompp -f ../md.mdp -c npt.gro -r npt.gro -t npt.cpt -p ../topol.top -o md.tpr
   gmx_gpu mdrun -v -deffnm md
fi

if [ -f md.xtc ]; then
    # Removing PBC effects, by placing the center of mass of the molecules in the box, centering the system using the protein, and placing all atoms at the closest distance from the center of the box (useful for visualization of truncated octahedron).  
    gmx_gpu trjconv -s md.tpr -f md.xtc -o md_center.xtc -pbc mol -center -ur compact <<< "1 0"
fi
