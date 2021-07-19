# This file is a part of the TRain program
# Author: Austin Seamann & Dario Ghersi
# Version: 1.0
# Last Updated: July 19th, 2021

import argparse
import subprocess
import time
import os
from shutil import copyfile


####################
# Global Variables #
####################

version = "linuxgccrelease"
program_dir = os.getcwd()
rosetta_dir = "pizza"


####################
#     Methods      #
####################
# Method: prep_numbers()
# Goal: Determine number of cpus needed for each run, rigid & flexible.
# Input:
#   cores - maximum number of cores allotted (both)
#   flex - number of flexible runs (rigid)
#   docking - number of docking runs (both)
#   refine - number of refinement runs (both)
#   pmhc - number of pmhc relax runs (flexible)
#   xml - number of xml relax runs (flexible)
#   bb - number of back bone rubs runs (flexible)
#   fast - number of fast relax runs (flexible)
# Output:
#   run_info: dictionary of parameters needed for running docking protocols
def prep_numbers(cores, flex, relax, docking, refine, pmhc, xml, bb, fast):
    run_info = {}
    # Determines docking core count
    run_info["docking"] = docking
    if cores > docking:
        run_info["cpu_docking"] = docking + 1
    else:
        run_info["cpu_docking"] = cores
    # Determines refinement core count
    run_info["refine"] = refine
    if cores > refine:
        run_info["cpu_refine"] = refine + 1
    else:
        run_info["cpu_refine"] = cores
    # Flexible docking numbers
    if flex:
        run_info["pmhc"] = pmhc  # Number of pmhc relax runs
        # Determine number of cpu cores for pmhc relax
        if cores > pmhc:
            run_info["cpu_pmhc"] = pmhc + 1  # Plus 1 for controller core
        else:
            run_info["cpu_pmhc"] = cores  # Use all cores
        run_info["xml"] = xml # number of xml relax runs
        # Determine number of cpu cores for xml relax
        if cores > xml:
            run_info["cpu_xml"] = xml + 1
        else:
            run_info["cpu_xml"] = cores
        run_info["bb"] = bb  # Number of back bone rub relax runs
        if cores > bb:
            run_info["cpu_bb"] = bb + 1
        else:
            run_info["cpu_bb"] = cores
        run_info["fast"] = fast  # Number of fast relax runs
        if cores > fast:
            run_info["cpu_fast"] = fast
        else:
            run_info["cpu_fast"] = cores
    # Rigid docking numbers
    else:
        run_info["relax"] = relax
        if cores > relax:
            run_info["cpu_relax"] = relax + 1
        else:
            run_info["cpu_relax"] = cores
    return run_info


# Method: run_rigid()
# Goal: Run rigid docking with either default or provided parameters from user.
# Input:
#   pdb: location of pdb(s)
#   single_file: boolean value of if a single pdb file or a directory of pdb files
def run_rigid(pdb, single_file):
    return "help"


# Method: run_flexible()
# Goal: Run flexib le docking with either default or provided parameters from user.
# Input:
#   pdb: location of pdb(s)
#   single_file: boolean value of if a single pdb file or a directory of pdb files
def run_flexible(pdb, single_file):
    return "help"


####################
#     Controls     #
####################
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="PDB file or folder of PDBs for docking", type=str)
    parser.add_argument("-i", "--initialize", help="Location of rosetta directory", type=str,
                        default="default")
    parser.add_argument("-l", "--linux", help="Changes to linux runnable program", default=True,
                        dest='linux', action='store_true')
    parser.add_argument("-m", "--mac", help="Changes to mac runnable program", default=False,
                        dest='mac', action='store_true')
    parser.add_argument("-r", "--rigid", help="Initializes flexible docking", default=True,
                        action="store_true")
    parser.add_argument("-f", "--flexible", help="Initializes flexible docking", default=False,
                        action='store_true')
    parser.add_argument("-c", "--cores", help="Max number of cpu cores provided, default=all",
                        default=os.cpu_count(), type=int)
    parser.add_argument("-a", "--relax", help="(Rigid) Number of relax runs performed", default=100,
                        type=int)
    parser.add_argument("-d", "--docking", help="Number of docking runs performed", default=10000,
                        type=int)
    parser.add_argument("-p", "--pmhc", help="(Flex) Number of pmhc relax runs", default=100,
                        type=int)
    parser.add_argument("-x", "--xml", help="(Flex) Number of xml relax runs for TCR", default=40,
                        type=int)
    parser.add_argument("-b", "--bb", help="(Flex) Number of back bone rub runs for TCR", default=30,
                        type=int)
    parser.add_argument("-s", "--fast", help="(Flex) Number of fast relax runs for TCR", default=30,
                        type=int)
    parser.add_argument("-e", "--refine", help="Number of refinement runs done, post dock", default=100,
                        type=int)
    return parser.parse_args()


def main():
    args = parse_args()
    # Initializing rosetta folder
    if args.initialize != "default":
        with open(__file__, "r") as f:
            lines = f.read().split('\n')
            with open("test.py", "w") as f1:
                for line in lines:
                    if line.startswith('rosetta_dir = "'):
                        line = line.split('"')[0] + '"' + args.initialize + '"'
                    f1.write(line + "\n")
    else:
        # Initialize version of Linux or Mac release
        global version
        if args.mac:
            version = "macosclangrelease"
        elif args.linux:
            version = "linuxgccrelease"
        # Initialize location of rosetta install
        global rosetta_dir
        rosetta_dir = args.rosetta
        # Initialize flexible or rigid docking & determine if a batch of pdbs or a single file
        run_info = prep_numbers(args.cores, args.flexible, args.relax, args.docking, args.refine, args.pmhc, args.xml,
                                args.bb, args.fast)
        if args.flexible:
            run_flexible(args.pdb, os.path.isfile(args.pdb), run_info)
        elif args.rigid:
            run_rigid(args.pdb, os.path.isfile(args.pdb), run_info)


if __name__ == '__main__':
    main()

