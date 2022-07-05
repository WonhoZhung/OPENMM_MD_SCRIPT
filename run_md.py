import os 
import sys
import openmm
from time import time
from rdkit import Chem
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import SystemGenerator, GAFFTemplateGenerator
from openmm import app, Platform, LangevinIntegrator
from openmm.app import PDBFile, Simulation, Modeller, PDBReporter, DCDReporter, StateDataReporter, ForceField
import md_params


def set_platform():
    speed = 0
    for i in range(Platform.getNumPlatforms()):
        p = Platform.getPlatform(i)
        p_speed = p.getSpeed()
        if p_speed > speed:
            platform = p
            speed = p_speed

    if platform.getName() == 'CUDA':
        platform.setPropertyDefaultValue('Precision', 'mixed')
    return platform

def prepare_system(ligand_fn, protein_fn, complex_fn=None, solvate=False, verbose=False):
    ligand_mol = Molecule(Chem.SDMolSupplier(ligand_fn)[0])
    protein_mol = PDBFile(protein_fn)
    
    modeller = Modeller(protein_mol.topology, protein_mol.positions)
    modeller.add(ligand_mol.to_topology().to_openmm(), ligand_mol.conformers[0])

    if complex_fn is not None:
        with open(complex_fn, 'w') as w:
            PDBFile.writeFile(modeller.topology, modeller.positions, w)

    if verbose: 
        print("System generation starts!")
    st = time()
    system_generator = SystemGenerator(
                forcefields=md_params.FORCEFIELDS,
                small_molecule_forcefield=md_params.SMALL_MOLECULE_FORCEFIELD,
                forcefield_kwargs=md_params.FORCEFIELD_KWARGS,
                molecules=[ligand_mol]
    )

    if solvate:
        modeller.addSolvent(
                    system_generator.forcefield,
                    model='tip3p',
                    padding=md_params.PADDING,
                    ionicStrength=md_params.IONIC_STRENGTH,
                    positiveIon=md_params.POSITIVE_ION,
                    negativeIon=md_params.NEGATIVE_ION
        )
    
        if complex_fn is not None:
            with open(complex_fn[:-4] + "_solv.pdb", 'w') as w:
                PDBFile.writeFile(modeller.topology, modeller.positions, w)

    system = system_generator.create_system(
                modeller.topology, 
                molecules=ligand_mol
    )
    integrator = LangevinIntegrator(
                md_params.TEMPERATURE,
                md_params.FRICTION_COEFF,
                md_params.STEP_SIZE
    )
    if solvate:
        system.addForce(
                    openmm.MonteCarloBarostat(
                            md_params.PRESSURE, 
                            md_params.TEMPERATURE, 
                            25
                    )
        )
    if verbose: 
        print(f"System generation ends! --- Duration: {time() - st:.2f}(s)")
    return modeller, system, integrator

def run_simulation(modeller, system, integrator, platform, traj_fn, log_fn, verbose=False):
    simulation = Simulation(
            topology=modeller.topology,
            system=system,
            integrator=integrator,
            platform=platform
    )

    # 1. Energy minimization
    if verbose: 
        print("Energy minimization starts!")
    st = time()
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    if verbose: 
        print(f"Energy minimization ends! --- Duration: {time() - st:.2f}(s)")

    # 2. Langevin NVT equilibrium
    if verbose: 
        print("Equilibrium starts!")
        st = time()
    simulation.context.setVelocitiesToTemperature(md_params.TEMPERATURE)
    simulation.step(md_params.NUM_EQ_STEPS)
    if verbose: 
        print(f"Equilibrium ends! --- Duration: {time() - st:.2f}(s)")

    # 3. MD simulation
    if verbose: 
        print("MD simulation starts!")
    st = time()
    #simulation.reporters.append(PDBReporter(traj_fn, \
    #        md_params.REPORTING_INTERVAL, enforcePeriodicBox=False))
    simulation.reporters.append(DCDReporter(traj_fn, \
            md_params.REPORTING_INTERVAL, enforcePeriodicBox=False))
    simulation.reporters.append(StateDataReporter(log_fn, \
            md_params.REPORTING_INTERVAL * 5, step=True, \
            potentialEnergy=True, temperature=True))
    simulation.step(md_params.NUM_MD_STEPS)
    if verbose: 
        print(f"MD simulation ends! --- Duration: {time() - st:.2f}(s)")
    return


