from pdbfixer import PDBFixer
from openmm.app import PDBFile
from glob import glob
from rdkit import Chem
import sys
import os
import md_params

def prepare_ligand(ligand_fn1, ligand_fn2):
    mol = Chem.SDMolSupplier(ligand_fn1)[0]
    for a in mol.GetAtoms():
        a.SetNumExplicitHs(a.GetNumRadicalElectrons())
        a.SetNumRadicalElectrons(0)
    mol_addh = Chem.AddHs(mol, addCoords=True)
    Chem.AssignAtomChiralTagsFromStructure(mol_addh)

    writer = Chem.SDWriter(ligand_fn2)
    writer.write(mol_addh)
    writer.close()
    return

def prepare_ligands(ligands_fn, ligand_dir, addH=True):
    mols = Chem.SDMolSupplier(ligands_fn)
    for mol in mols:
        if addH:
            for a in mol.GetAtoms():
                a.SetNumExplicitHs(a.GetNumRadicalElectrons())
                a.SetNumRadicalElectrons(0)
            mol_addh = Chem.AddHs(mol, addCoords=True)
            Chem.AssignAtomChiralTagsFromStructure(mol_addh)
        else:
            mol_addh = mol

        writer = Chem.SDWriter(ligand_dir + mol.GetProp("_Name") + '.sdf')
        writer.write(mol_addh)
        writer.close()
    return

def prepare_protein(protein_fn1, protein_fn2, pH=md_params.PH):
    fixer = PDBFixer(filename=protein_fn1)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH)
    
    fixer.removeHeterogens(keepWater=False)

    with open(protein_fn2, 'w') as w:
        PDBFile.writeFile(fixer.topology, fixer.positions, file=w, keepIds=True)
    return

# def run(key):
#     protein_fn1 = f""
#     protein_fn2 = f""
#     ligand_fn1 = f""
#     ligand_fn2 = f""
# 
#     prepare_protein(protein_fn1, protein_fn2)
#     prepare_ligand(ligand_fn1, ligand_fn2)
#     return

if __name__ == "__main__":

    import argparse
    from time import time
    
    parser = argparse.ArgumentParser()

    parser.add_argument("--p_in", help="protein input; PDB format", type=str, required=True)
    parser.add_argument("--p_out", help="protein output; PDB format", type=str, required=True)
    parser.add_argument("--l_in", help="ligand input; SDF format", type=str, required=True)
    parser.add_argument("--l_out", help="ligand output; SDF format", type=str, required=True)
    #parser.add_argument("-v", help="verbose", action="store_true")

    args = parser.parse_args()

    prepare_protein(args.p_in, args.p_out)
    prepare_ligand(args.l_in, args.l_out)
