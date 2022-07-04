from simtk import unit
from openmm import app


PH = 7.4
TEMPERATURE = 303.15 * unit.kelvin
PRESSURE = 1.0 * unit.atmosphere
IONIC_STRENGTH = 0.1 * unit.molar
POSITIVE_ION = "K+"
NEGATIVE_ION = "Cl-"
PADDING = 10.0 * unit.angstrom
FRICTION_COEFF = 1 / unit.picosecond
STEP_SIZE = 0.002 * unit.picoseconds # 2 fs
NUM_EQ_STEPS = int(5e4) # 100 ps
REPORTING_INTERVAL = int(5e3) # 10 ps
NUM_MD_STEPS = int(5e5) # 1 ns

FORCEFIELD_KWARGS = {
            "constraints": app.HBonds,
            "rigidWater": True,
            "removeCMMotion": False,
            "hydrogenMass": 4 * unit.amu
}
FORCEFIELDS = ['amber/ff14SB.xml', 'amber/tip3p_standard.xml']
SMALL_MOLECULE_FORCEFIELD = 'gaff-2.11'
