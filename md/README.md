# README - Molecular Dynamics Files

Directory contents:

- `inputs`: GROMACS simulation inputs for all evaluated systems (`gro`, `top`, `mdp`), generated as follows:

  ```bash
  # Obtaining \lstinline|gro| from \lstinline|pdb| file, output from AlphaFold2:
  gmx pdb2gmx -f protein.pdb -o protein.gro -ter -ignh
  # Select force field, water model and terminal types (zwitterionic)
  # Force field parameters must be present in the directory, see README.md in the project
  # root for respective links where they can be obtained

  # Preparing the simulation box:
  # Truncated octahedron box with 1.5 nm padding
  gmx editconf -f protein.gro -o protein_box.gro -c -d 1.5 -bt octahedron

  # Solvation with CHARMM36m:
  gmx solvate -cp protein_box.gro -cs spc216.gro -o protein_water.gro -p topol.top

  # OR:

  # Solvation with a99SB-disp:
  gmx solvate -cp protein_box.gro -cs ./a99SBdisp.ff/a99SBdisp_water.gro -o protein_water.gro -p topol.top

  # Neutralization and correcting ionic strength:
  gmx grompp -f ions.mdp -c protein_water.gro -p topol.top -o ions.tpr
  # Replaces water molecules with Na+ and Cl- setting concentration to 0.15 M (for neutralization only, omit "-conc 0.15")
  gmx genion -s ions.tpr -o protein_neutral.gro -p topol.top -pname NA -nname CL -conc 0.15 -neutral <<< "13"
  ```

- `bash_scripts`: `.sh` scripts used to lauch simulations (from files in `inputs`), and to extend simulations up to 5 μs
- `data`: processed simulation timeseries for all evaluated systems. Radius of gyration as a function of time is found in `rg.csv` files, secondary structure content estimated with the DSSP algorithm in `dssp.npy` files, and chemical shifts predicted with SPARTA+ in `shifts_merged.pkl` files.
- `analysis`: various Python scripts used to analyze MD trajectories
  - `props.py`: obtaining `rg.csv`, `dssp.npy`, chemical shifts (in chunks) from `.xtc` files (requires MDTraj and SPARTA+)
  - `utils.py`: utility functions for formatting chemical shift files and merging them to obtain `shifts_merged.pkl`
  - `io.py`: functions for loading simulation data and formatting into structured `pandas.DataFrame` objects (requires pandas and numpy packages)
  - `equilibration_time.py`: contains functions for automated determination of a trajectory's equilibration time based on estimated relaxation times of various properties
  - `blocking.py`: contains functions used for calculating aggregate statistics and standard errors for simulation-derived properties, while accounting for their autocorrelation using blocking analysis
