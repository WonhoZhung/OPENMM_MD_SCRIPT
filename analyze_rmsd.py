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
)

atoms = t.topology.select("chainid 1")
rmsds_lig = md.rmsd(t, t, frame=0, atom_indices=atoms, parallel=True, precentered=False)
#atoms = t.topology.select("chainid 0 and backbone")
#rmsds_bck = md.rmsd(t, t, frame=0, atom_indices=atoms, parallel=True, precentered=False)

#df = pd.DataFrame({"time": t.time, "lig": np.array(rmsds_lig), "bck": np.array(rmsds_bck)})
df = pd.DataFrame({"time": t.time * 10, "lig": np.array(rmsds_lig)})

sns.set_theme(style="darkgrid")
sns.lineplot(data=df, x="time", y="lig", label="Ligand")
#sns.lineplot(data=df, x="time", y="bck", label="Backbone")

plt.xlabel('Time (ps)')
plt.ylabel('RMSD (Å)')

plt.show()

print(f"Mean RMSD (Å): {np.array(rmsds_lig).mean():.3f}")
print(f"Std RMSD (Å): {np.array(rmsds_lig).std():.3f}")
