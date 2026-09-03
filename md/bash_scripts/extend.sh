#!/bin/bash
OUTDIR="out"
cd $OUTDIR

if [ -f md_center.xtc ]; then
   rm md_center.xtc
fi

# extend until 5 microseconds
gmx_gpu convert-tpr -s md_ext.tpr -until 5000000 -o md_ext.tpr
gmx_gpu mdrun -s md_ext.tpr -cpi md.cpt  -v -deffnm md

if [ -f md.xtc ]; then
    # Re-center
    gmx_gpu trjconv -s md_ext.tpr -f md.xtc -o md_center.xtc -pbc mol -center -ur compact <<< "1 0"
    # Protein only traj
    gmx_gpu trjconv -s md_ext.tpr -f md_center.xtc -o md_protein.xtc <<< "1"
fi
