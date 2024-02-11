#!/usr/bin/python3

######################################################################
# TurnTable.py -- A component of TRain                               #
# Copyright: Austin Seamann & Dario Ghersi                           #
# Version: 1.0                                                       #
# Last Updated: January 12th, 2024                                   #
# Goal: Turn table takes the input of modeled or crystal structure   #
#       PDBs to construct of PDB with a TCR paired with a pMHC in a  #
#       specified manor for RosettaDock4.0                           #
#                                                                    #
# Named arguments: -t --tcr (TCR PDB file or folder of TCR PDBs)     #
#                  -p --pmhc (pMHC PDB file or folder of pMHC PDBs)  #
#                  -d --pdbs (Submit folder of PDBs and all TCRs and #
#                             will be swapped)                       #
#                  -a --trimA (Trim TCR alpha, needed if full        #
#                             crystal structure: variable + constant)#
#                  -b --trimB (Trim TCR beta, needed if full crystal #
#                             structure: variable + constant)        #
#                  -m --trimM (Prevent trimming of MHC)              #
######################################################################


import argparse
from util.PDB_Tools_V3 import PdbTools3
import os
from shutil import copyfile


####################
# Global Variables #
####################
tool = PdbTools3()  # Initialize PDB tools


#################
#    Methods    #
#################
def align_chains(tcr_file, pmhc_file, final_name):
    """
    Prepares cmd file and submits to run_chimera for aligning pdb to reference
    and producing a separated protein and peptide pdb.

    Parameters
    __________
    tcr_file : str
        location of tcr pdb
    pmhc_file : str
        location of pmhc pdb
    final_name : str
        name of resulting paired TCRpMHC file
    """
    reference_pdb = os.path.join(os.path.dirname(__file__), "reference.pdb")
    try:
        tool.set_file_name(reference_pdb)
    except:
        print("Missing reference pdb - backup file in /TRain/data/")
    tool.superimpose(tcr_file, "DE", "DE", "reference_tmp.pdb")
    tool.set_file_name(pmhc_file)
    tool.superimpose("reference_tmp.pdb", "A", "A", pmhc_file)
    tool.join(pmhc_file, tcr_file, final_name)
    os.remove("reference_tmp.pdb")


def tcr_pmhc_pair(tcr_dir, pmhc_dir, tcr_multi, pmhc_multi, trim_a, trim_b, trim_m):
    """
    Goals: Pair tcr(s) and pmhc(s) submitted. Determines if single file is submitted of either or directory of either.

    Parameters
    __________
    tcr_dic : str
        directory of the location of the TCR pdb or TCR pdb folder
    pmhc_dir : str
        directory of the location of the pMHC pdb or pMHC pdb folder
    tcr_multi : boolean
        if true, a folder of TCR pdbs has been submitted
    pmhc_multi : boolean
        if true, a folder of pMHC pdbs has been submitted
    trim_a : boolean
        if true, don't trim alpha chain to just the variable region
    trim_b : boolean
        if true, don't trim beta chain to just the variable region
    trim_m : boolean
        if true, don't trim the MHC chain
    """
    directory = os.getcwd()
    # Full directory locations of TCRs and pMHCs
    locations = {"TCR": [], "pMHC": []}
    if tcr_multi:
        for tcr in sorted(os.listdir(tcr_dir)):
            if tcr.endswith(".pdb"):
                locations["TCR"].append(os.getcwd() + "/" + tcr_dir + "/" + tcr)
    else:
        locations["TCR"].append(tcr_dir)
    if pmhc_multi:
        for pmhc in sorted(os.listdir(pmhc_dir)):
            if pmhc.endswith(".pdb"):
                locations["pMHC"].append(os.getcwd() + "/" + pmhc_dir + "/" + pmhc)
    else:
        locations["pMHC"].append(pmhc_dir)
    if "Paired" not in os.listdir():
        os.mkdir("Paired")
    for tcr_location in locations["TCR"]:
        if len(locations["TCR"]) <= 1:
            tmp_tcr = os.getcwd() + "/" + tcr_location.split("/")[-1].split(".")[0] + "_tmpt.pdb"
        else:
            tmp_tcr = "/".join(tcr_location.split("/")[:-1]) + "/" + tcr_location.split("/")[-1].split(".")[0] \
                           + "_tmpt.pdb"
        copyfile(tcr_location, tmp_tcr)
        tool.set_file_name(tmp_tcr)
        num_of_chains = len(tool.get_chains())
        if num_of_chains > 2:  # Assume crystal structure of TCR
            tool.clean_pdb()  # Updates to A, B, C, D, E format of PDB and removes any additional chains
            tool.split_tcr(tmp_tcr)  # Splits to a tmp file for just a/b chains
        else:
            chain_order = tool.get_chains()
            if chain_order == ["B", "A"]:  # RepBuilder output correction
                tool.reorder_chains(["A", "B"])
            tool.update_label({'A': 'D', 'B': 'E'})
            tool.renumber_docking()
        if trim_a:  # Trimming of alpha chain
            tool.trim_chain("D", 107)  # Trim Alpha at 107
        else:
            if len(tool.get_amino_acid_on_chain("D")) >= 120:
                print("Trim A flag not called, suggested for this structure")
        if trim_b:  # Trimming of beta chain
            tool.trim_chain("E", 113)  # Trim Beta at 113
        else:
            if len(tool.get_amino_acid_on_chain("E")) >= 120:
                print("Trim B flag not called, suggested for this structure")
        tool.center(tmp_tcr)
        for pmhc_location in locations["pMHC"]:
            if len(locations["pMHC"]) <= 1:
                tmp_pmhc = os.getcwd() + "/" + pmhc_location.split("/")[-1].split(".")[0] + "_tmpp.pdb"
            else:
                tmp_pmhc = "/".join(pmhc_location.split("/")[:-1]) + "/" + pmhc_location.split("/")[-1].split(".")[0] \
                           + "_tmpp.pdb"
            copyfile(pmhc_location, tmp_pmhc)
            tool.set_file_name(tmp_pmhc)
            num_chains = len(tool.get_chains())
            if num_chains > 3:  # assume it's just the pMHC component
                tool.clean_pdb()  # Updates to A, B, C, D, E format of PDB and removes any additional chains
            tool.split_pmhc(tmp_pmhc)  # Splits to a tmp file for peptide and mhc
            if not trim_m:  # Prevent trimming of MHC chain
                tool.trim_chain("A", 181)  # Trim MHC at 181
            # Name of new PDB
            new_name = "_".join(tmp_tcr.split("/")[-1].split("_")[:-1]) + "_" \
                       + "_".join(tmp_pmhc.split("/")[-1].split("_")[:-1]) + ".pdb"
            name1 = directory + "/Paired/" + new_name  # TCR first, pMHC second
            align_chains(tmp_tcr, tmp_pmhc, name1)  # Pair TCRs to pMHCs
            tool.set_file_name(directory + '/Paired/' + new_name)
            tool.renumber_docking()  # Renumber after creating paired TCRpMHCs
            os.remove(tmp_pmhc)
        os.remove(tmp_tcr)


####################
#     Controls     #
####################
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tcr", help="TCR pdb file or folder of TCR PDBs", type=str, default="...")
    parser.add_argument("-p", "--pmhc", help="pMHC pdb file or folder of pMHC PDBs", type=str, default="...")
    parser.add_argument("-d", "--pdbs", help="Submit folder of pdbs and all TCRs and pMHCs will be swapped", type=str,
                        default="...")
    parser.add_argument("-a", "--trimA", help="Prevent trimming TCR alpha, needed if full crystal structure",
                        action="store_true", default=False)
    parser.add_argument("-b", "--trimB", help="Prevent trimming TCR beta, needed if full crystal structure",
                        action="store_true", default=False)
    parser.add_argument("-m", "--trimM", help="Prevent trimming of MHC", action="store_true", default=False)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.tcr != "...":
        tcr_multi = False  # True if directory of tcrs and not single file
        pmhc_multi = False  # True if directory of pmhcs and not single file
        if not os.path.isfile(args.tcr):
            tcr_multi = True
        tcr_location = args.tcr
        if args.pmhc != "...":
            if not os.path.isfile(args.pmhc):
                pmhc_multi = True
            pmhc_location = args.pmhc
            tcr_pmhc_pair(tcr_location, pmhc_location, tcr_multi, pmhc_multi, args.trimA, args.trimB, args.trimM)
        else:
            print("Please provide pMHC files.")
    # If directory of several PDBs and to swap out each antigen and each TCR
    elif args.pdbs != "...":
        tcr_pmhc_pair(args.pdbs, args.pdbs, True, True, args.trimA, args.trimB, args.trimM)


if __name__ == '__main__':
    main()
