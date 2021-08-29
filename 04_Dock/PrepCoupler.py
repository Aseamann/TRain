# This file is a part of the TRain program
# Author: Austin Seamann & Dario Ghersi
# Version: 1.0
# Last Updated: August 29th, 2021

from PDB_Tools_V3 import PdbTools3
import os

####################
# Global Variables #
####################
tool = PdbTools3()  # Initialize PDB tools

#################
#    Methods    #
#################


####################
#     Controls     #
####################
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="PDB file or folder of PDBs for docking", type=str)
    return parser.parse_args()


def main():
    args = parse_args()
    if os.path.isfile(args.pdb):
        for pdb in sorted(os.listdir(os.getcwd)):
            if pdb.endswith(".pdb"):
                tool.set_file_name(pdb)
                tool.clean_pdb()  # Updates to A, B, C, D, E format of PDB and removes any additional chains
                tool.renumber_docking()  # Update numbering of PDB to Rosetta standard
    else:
        tool.set_file_name(args.pdb)
        tool.clean_pdb()  # Updates to A, B, C, D, E format of PDB and removes any additional chains
        tool.renumber_docking()  # Update numbering of PDB to Rosetta standard


if __name__ == "__main__":
    main()
