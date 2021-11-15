# This file is a part of the TRain program
# Author: Austin Seamann & Dario Ghersi
# Version: 0.1
# Last Updated: November 14th, 2021

import argparse
from PDB_Tools_V3 import PdbTools3
import os


####################
#     Controls     #
####################
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="PDB file or folder of PDBs for docking", type=str)
    return parser.parse_args()


def main():
    args = parse_args()
    tool = PdbTools3()  # Initialize PDB tools
    if os.path.isdir(args.pdb):
        for pdb in sorted(os.listdir(os.getcwd() + "/" + args.pdb)):
            if pdb.endswith(".pdb"):
                file_name = os.getcwd() + "/" + args.pdb + "/" + pdb
                tool.set_file_name(file_name)
                tool.clean_pdb()  # Updates to A, B, C, D, E format of PDB and removes any additional chains
                print(tool.get_chains())
                tool.remove_chain("B")  # Remove B2M
                tool.trim_chain("D", 107)  # Trim Alpha at 107
                tool.trim_chain("E", 113)  # Trim Beta at 113
                tool.trim_chain("A", 181)  # Trim MHC at 181
                tmp_tcr = ".".join(file_name.split(".")[:-1]) + "_tmpt.pdb"
                tool.split_tcr(tmp_tcr, True)
                tool.set_file_name(tmp_tcr)
                tool.center(tmp_tcr)
                tool.set_file_name(file_name)
                tool.superimpose(tmp_tcr, "DE", "DE", file_name)
                os.remove(tmp_tcr)
                tool.renumber_docking()  # Update numbering of PDB to Rosetta standard
    else:
        tool.set_file_name(args.pdb)
        tool.clean_pdb()  # Updates to A, B, C, D, E format of PDB and removes any additional chains
        tool.remove_chain("B")  # Remove B2M
        tool.trim_chain("D", 107)  # Trim Alpha at 107
        tool.trim_chain("E", 113)  # Trim Beta at 113
        tool.trim_chain("A", 181)  # Trim MHC at 181
        tmp_tcr = ".".join(args.pdb.split(".")[:-1]) + "_tmpt.pdb"
        tool.split_tcr(tmp_tcr, True)
        tool.set_file_name(tmp_tcr)
        tool.center(tmp_tcr)
        tool.set_file_name(args.pdb)
        tool.superimpose(tmp_tcr, "DE", "DE", args.pdb)
        os.remove(tmp_tcr)
        tool.renumber_docking()  # Update numbering of PDB to Rosetta standard


if __name__ == "__main__":
    main()
