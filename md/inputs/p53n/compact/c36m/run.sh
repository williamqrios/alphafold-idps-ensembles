#!/bin/bash
REPLICA="long"
cd $REPLICA

gmx grompp -f ../em.mdp -c ../protein_neutral.gro -p ../topol.top -o em.tpr
gmx mdrun -v -deffnm em

# Equilibration NVT for 1 ns (full restraint)
gmx grompp -f ../nvt.mdp -c em.gro -r em.gro -p ../topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt

# Equilibration NPT for 1 ns (full restraint), -maxwarn required due to Berendsen barostat 
if [ -f nvt.xtc ]; then
    gmx grompp -f ../npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p ../topol.top -o npt.tpr -maxwarn 2
    gmx mdrun -v -deffnm npt
fi

# Production
if [ -f npt.xtc ]; then
   gmx grompp -f ../md.mdp -c npt.gro -r npt.gro -t npt.cpt -p ../topol.top -o md.tpr
   gmx mdrun -v -deffnm md
fi

if [ -f md.xtc ]; then
    # Removing PBC effects, by placing the center of mass of the molecules in the box, centering the system using the protein, and placing all atoms at the closest distance from the center of the box (useful for visualization of truncated octahedron).  
    gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -pbc mol -center -ur compact <<< "1 0"
    #gmx trjconv -s md.tpr -f md_center.xtc -o at_0ns.gro -pbc mol -center -dump 0 <<< "1 0"
    #gmx grompp -f md.mdp -c at_0ns.gro -p topol.top -o at_0ns.tpr
    #gmx gyrate -s at_0ns.tpr -f md_center.xtc -o gyrate_0_200ns.xvg -tu ns <<< "1"
fi
