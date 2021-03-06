import mdtraj as md
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import sys
import argparse
from time import time

parser = argparse.ArgumentParser()

parser.add_argument("-t", help="trajectory file; DCD format", type=str, required=True)
parser.add_argument("-c", help="system topology file; PDB format", type=str, required=True)

args = parser.parse_args()

t = md.load(
    args.t,
    top=args.c
) # trajectory
r = md.load(
    args.c
) # reference

t.image_molecules()
t.save(args.t[:-4] + "_center.dcd")

lig_atoms = t.topology.select("chainid 1")
bck_atoms = t.topology.select("chainid 0 and backbone")
t.superpose(r, atom_indices=bck_atoms)
rmsds_lig = np.sqrt(3*np.mean((t.xyz[:, lig_atoms, :] - t.xyz[0:1, lig_atoms, :])**2, axis=(1,2)))

df = pd.DataFrame({"time": t.time, "lig": np.array(rmsds_lig)})

sns.set_theme(style="darkgrid")
sns.lineplot(data=df, x="time", y="lig", label="Ligand")

plt.xlabel('Time steps')
plt.ylabel('RMSD (Å)')

plt.show()

print(f"Mean RMSD (Å): {np.array(rmsds_lig).mean():.3f}")
print(f"Std RMSD (Å): {np.array(rmsds_lig).std():.3f}")
