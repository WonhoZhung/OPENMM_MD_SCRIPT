import os
import sys
import argparse
from run_md import set_platform, prepare_system, run_simulation


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-p", help="protein_fn in PDB format", type=str)
    parser.add_argument("-l", help="ligand_fn in SDF format", type=str)
    parser.add_argument("-c", help="complex_fn in PDB format", type=str)
    parser.add_argument("-t", help="traj_fn in DCD format", type=str)
    parser.add_argument("-s", help="solvate system with water", action="store_true")
    parser.add_argument("-o", help="log_fn", type=str)
    parser.add_argument("-v", help="verbose", action="store_true")

    args = parser.parse_args()

    platform = set_platform()
    modeller, system, integrator = \
            prepare_system(args.l, args.p, args.c, args.s, args.v) 
    run_simulation(modeller, system, integrator, platform, args.t, args.o, args.v)
    return


if __name__ == "__main__":

    main()
