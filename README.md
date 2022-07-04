# OPENMM_MD_SCRIPT

## 1. Preparing input files
```
usage: prepare_input.py [-h] --p_in P_IN --p_out P_OUT --l_in L_IN --l_out L_OUT

optional arguments:
  -h, --help     show this help message and exit
  --p_in P_IN    protein input; PDB format
  --p_out P_OUT  protein output; PDB format
  --l_in L_IN    ligand input; SDF format
  --l_out L_OUT  ligand output; SDF format
```

## 2. Run MD
```
usage: main.py [-h] [-p P] [-l L] [-c C] [-t T] [-s] [-o O] [-v]

optional arguments:
  -h, --help  show this help message and exit
  -p P        protein_fn in PDB format
  -l L        ligand_fn in SDF format
  -c C        complex_fn in PDB format
  -t T        traj_fn in DCD format
  -s          solvate system with water
  -o O        log_fn
  -v          verbose
```

## 3. Analyze RMSD
```
usage: analyze_rmsd.py [-h] -t T -c C

optional arguments:
  -h, --help  show this help message and exit
  -t T        trajectory file; DCD format
  -c C        system topology file; PDB format
```
