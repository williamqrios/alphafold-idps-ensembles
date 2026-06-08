# AlphaFold2-generated ensembles of intrinsically disordered proteins (IDPs)
Directory structure:
```
ensembles/
├─ asyn/
│  ├─ deep_recycle/
│  ├─ shallow_dropout/
├─ hst5/
│  ├─ deep_recycle/
│  ├─ shallow_dropout/
├─ p53n/
│  ├─ deep_recycle/
│  ├─ shallow_dropout/
├─ prta/
│  ├─ deep_recycle/
│  ├─ shallow_dropout/
rgs/
selected/
scripts/
```
- `ensembles/`: directory containing 200 AlphaFold2-generated structures (`pdb` files) for the listed IDPs: alpha-synuclein (asyn), histatin 5 (hst5), p53-NTD (p53n), prothymosin alpha (prta). The 200 structures are categorized into two directories labeled according to AlphaFold2 inference parameters (deep recycle, shallow dropout).
- `rgs/`: directory containing calculated radii of gyration of generated structures as `csv` files (one `csv` file per protein/parameter combination).  
- `selected/`: directory containing 3 `pdb` files per protein, selected from the 200 structural ensemble, according to radius of gyration (Rg). Compact = lowest Rg, extended = highest Rg, exp = nearest to reported experimental value. These structures were used directly as starting conformations for all-atom molecular dynamics simulations.
- `scripts/`: directory containing bash (`sh`) scripts for generating AlphaFold structural ensembles using LocalColabFold, accompanied by the corresponding `fasta` files. See the [https://github.com/YoshitakaMo/localcolabfold](LocalColabFold repository) for more details on how to install LocalColabFold.