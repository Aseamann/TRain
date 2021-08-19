# This file is a part of the TRain program
# Author: Austin Seamann & Dario Ghersi
# Version: 0.2
# Last Updated: Aug 4th, 2021
import argparse
from PDB_Tools_V3 import PdbTools3
import os
import subprocess


####################
# Global Variables #
####################
chimera_install = "/Applications/Chimera.app/Contents/MacOS/chimera"  # chimera installation location
tool = PdbTools3()  # Initialize PDB tools


#################
# Methods
#################
# Method: run_chimera()
# Goal: Runs chimera commands based on submitted cmd file and list of pdbs to open
def run_chimera(cmd_file, open_pdb=[]):
    global chimera_install
    chimera_dir = chimera_install  # Default location of chimera on mac
    run_list = [chimera_dir, "--nogui", "--silent"]  # Needed base arguments
    for each in open_pdb:  # Fill run_list with list of PDBs to open
        run_list.append(each)
    run_list.append(cmd_file)  # Adds in the cmd_file
    subprocess.run(run_list)  # Runs chimera


# Method: align_chains()
# Goal: Prepares cmd file and submits to run_chimera for aligning pdb to reference
#       and producing a separated protein and peptide pdb.
def align_chains(tcr_file, pmhc_file, final_name):
    with open("prepare.cmd", "w") as cmd:
        cmd.write("mmaker #0.1 #1\n")
        cmd.write("close #0.1 #0.2\n")
        cmd.write("mmaker #0.3 #2\n")
        cmd.write("close #0.3 #0.4\n")
        cmd.write("combine #\n")
        cmd.write("close #1 #2\n")
        cmd.write("write #0 " + final_name + "\n")
        cmd.write("close #\n")
    pdb_list = ["reference.pdb", pmhc_file, tcr_file]
    run_chimera("prepare.cmd", pdb_list)


# Method: multiple_pdbs()
# Goals: If collection of TCRpMHC files, swap each TCR with each pMHC. Along with trim MHC, Alpha, and Beta then remove
#        b2m.
def multiple_pdbs(pdb_folder):
    directory = os.getcwd()  # Collect cwd
    # Loops through all PDBs in full directory
    # Prepares full TCR/pMHC files.
    for pdb in sorted(os.listdir(pdb_folder)):
        if pdb.endswith(".pdb"):
            tool.set_file_name(pdb_folder + pdb)  # Sets PDB in tool
            print(pdb)
            # TODO: Update this method to relabel based on input
            tool.clean_pdb()  # Updates to A, B, C, D, E format of PDB and removes any additional chains
            tool.remove_chain(tool.get_b2m_chain())  # Remove B2M
            tool.trim_chain("A", 181)  # Trim MHC at 181
            tool.trim_chain("D", 107)  # Trim Alpha at 107
            tool.trim_chain("E", 113)  # Trim Beta at 113
            tool.split_chains("DE", "_tcr", directory + '/TCRs/')  # Splits TCR
            tool.split_chains("AC", "_pmhc", directory + '/pMHCs/')  # Splits pMHC
    for pdb1 in sorted(os.listdir(pdb_folder)):
        if pdb1.endswith(".pdb"):
            for pdb2 in sorted(os.listdir(pdb_folder)):
                if pdb2.endswith(".pdb"):
                    tcr1 = directory + "/TCRs/" + pdb1[:4] + "_tcr.pdb"
                    pmhc2 = directory + "/pMHCs/" + pdb2[:4] + "_pmhc.pdb"
                    new_name = pdb1[:4] + "_" + pdb2[:4] + ".pdb"  # Name of new PDB
                    name1 = directory + "/Paired/" + new_name  # TCR first, pMHC second
                    tool.set_file_name(directory + '/Paired/' + new_name)
                    align_chains(tcr1, pmhc2, name1)  # Pair TCRs to pMHCs, pMHCs get opened first
                    tool.renumber_docking()  # Renumber after creating paired TCRpMHCs


# Method: tcr_pmhc_pair()
# Goals: Pair tcr(s) and pmhc(s) submitted. Determines if single file is submitted of either or directory of either.
def tcr_pmhc_pair(tcr_dir, pmhc_dir, tcr_multi, pmhc_multi):
    directory = os.getcwd()
    # Full directory locations of TCRs and pMHCs
    locations = {"TCR": [], "pMHC": []}
    if tcr_multi:
        for tcr in sorted(os.listdir(tcr_dir)):
            if tcr.endswith(".pdb"):
                locations["TCR"].append(tcr)
    else:
        locations["TCR"].append(tcr_dir)
    if pmhc_multi:
        for pmhc in sorted(os.listdir(pmhc_dir)):
            if pmhc.endswith(".pdb"):
                locations["pMHC"].append(pmhc)
    else:
        locations["pMHC"].append(pmhc_dir)
    os.mkdir("Paired")
    for tcr_location in locations["TCR"]:
        tool.set_file_name(tcr_location)
        tool.clean_pdb()  # Updates to A, B, C, D, E format of PDB and removes any additional chains
        tmp_tcr = tcr_location.split("/")[-1].split(".")[0] + ".tmp"
        tool.split_tcr(tmp_tcr)  # Splits to a tmp file for just a/b chains
        for pmhc_location in locations["pMHC"]:
            tool.set_file_name(pmhc_location)
            tool.clean_pdb()  # Updates to A, B, C, D, E format of PDB and removes any additional chains
            tmp_pmhc = pmhc_location.split("/")[-1].split(".")[0] + ".tmp"
            tool.split_pmhc(tmp_pmhc)  # Splits to a tmp file for peptide and mhc
            new_name = tmp_tcr.split(".")[0] + "_" + tmp_pmhc.split(".")[0] + ".pdb"  # Name of new PDB
            name1 = directory + "/Paired/" + new_name  # TCR first, pMHC second
            tool.set_file_name(directory + '/Paired/' + new_name)
            align_chains(tmp_tcr, tmp_pmhc, name1)  # Pair TCRs to pMHCs, pMHCs get opened first
            tool.renumber_docking()  # Renumber after creating paired TCRpMHCs
            os.remove(tmp_pmhc)
        os.remove(tmp_tcr)


####################
#     Controls     #
####################
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tcr", help="TCR pdb file or folder of TCR PDBs", type=str, default="...")
    parser.add_argument("-m", "--pmhc", help="pMHC pdb file or folder of pMHC PDBs", type=str, default="...")
    parser.add_argument("-p", "--pdbs", help="Submit folder of pdbs and all TCRs and pMHCs will be swapped", type=str,
                        default="...")
    parser.add_argument("-c", "--chimera", help="(Required) Chimera install location, required before first run.", type=str,
                        default="/Applications/Chimera.app/Contents/MacOS/chimera")
    parser.add_argument("-l", "--label", help="Labeling update ex. ABCDE D=Alpha E=Beta", type=str, default="ABCDE")
    parser.add_argument("--b2m", help="Prevent removal of B2M", action="store_true", default=False)
    parser.add_argument("--trimA", help="Prevent trimming of TCR alpha", action="store_true", default=False)
    parser.add_argument("--trimB", help="Prevent trimming of TCR beta", action="store_true", default=False)
    parser.add_argument("--trimM", help="Prevent trimming of MHC", action="store_true", default=False)
    return parser.parse_args()


def main():
    args = parse_args()
    # Initializing chimera
    if args.chimera != "/Applications/Chimera.app/Contents/MacOS/chimera":
        with open(__file__, "r") as f:
            lines = f.read().split('\n')
            with open(__file__, "w") as f1:
                for line in lines:
                    if line.startswith('chimera_install = "'):
                        line = line.split('"')[0] + '"' + args.initialize + '"'
                    f1.write(line + "\n")
    if args.tcr != "...":
        tcr_multi = False  # True if directory of tcrs and not single file
        pmhc_multi = False  # True if directory of pmhcs and not single file
        if args.tcr.path.isfile():
            tcr_multi = True
        tcr_location = args.tcr
        if args.pmhc != "...":
            if args.pmhc.path.istfile():
                pmhc_multi = True
            pmhc_location = args.pmhc
            tcr_pmhc_pair(tcr_location, pmhc_location, tcr_multi, pmhc_multi)
        else:
            print("Please provide pMHC files.")
    # If directory of several PDBs and to swap out each antigen and each TCR
    if args.pdbs != "...":
        pdb_folder = args.pdbs  # Default structures
        multiple_pdbs(pdb_folder)


if __name__ == '__main__':
    main()
