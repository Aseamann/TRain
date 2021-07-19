import subprocess
import sys
import os
import argparse
from PDB_Tools_V3 import PdbTools3


verbose = False


# Manages program
def run_interface(tcr_dir, program, chains):
    global verbose
    tool = PdbTools3()  # Initialize PDB tools
    write_flag(chains)  # Creates Alpha and Beta flag files
    results = {}  # Results dictionary that will be used to create output ex. PDBid: {alpha: [scores], beta: [scores]}
    # Score = dG_separated
    # Loops through each PDB in submitted directory
    for pdb in os.listdir(tcr_dir):
        if pdb.endswith(".pdb"):
            tool.set_file_name(tcr_dir + "/" + pdb)
            results[pdb] = {"ALPHA": -1.0, "BETA": -1.0}
            tool.mute_aa(0, 1000, "B")  # muting the beta chain
            # Alpha
            if verbose:
                subprocess.run([program, "-s", tcr_dir + "/" + pdb, "@ALPHA_flags"])  # Runs InterfaceAnalyzer
            else:
                subprocess.run([program, "-s", tcr_dir + "/" + pdb, "@ALPHA_flags"],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # Runs InterfaceAnalyzer
            tool.unmute_aa(0, 1000, "B")  # un-muting beta chain
            tool.mute_aa(0, 1000, "A")  # muting the alpha chain
            # Beta
            if verbose:
                subprocess.run([program, "-s", tcr_dir + "/" + pdb, "@BETA_flags"])  # Runs InterfaceAnalyzer
            else:
                subprocess.run([program, "-s", tcr_dir + "/" + pdb, "@BETA_flags"],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # Runs InterfaceAnalyzer
            tool.unmute_aa(0, 1000, "A")  # un-muting alpha chain
            for score_file in os.listdir(os.getcwd()):  # Collecting scores
                if score_file.endswith(".sc"):
                    # Input to collect scores from output files
                    express = "tr -s ' ' < " + score_file + " | tail -1 | cut -f6 -d' ' | tr ' ' '\t'"
                    if score_file == "tcr_ALPHA.sc":
                        temp_score = subprocess.run(express, shell=True, stdout=subprocess.PIPE)  # Runs expression
                        results[pdb]["ALPHA"] = float(temp_score.stdout.decode('utf-8')[1:])  # Saves score to results
                        os.remove(score_file)  # Clears file before next run
                    if score_file == "tcr_BETA.sc":
                        temp_score = subprocess.run(express, shell=True, stdout=subprocess.PIPE)  # Runs expression
                        results[pdb]["BETA"] = float(temp_score.stdout.decode('utf-8')[1:])
                        os.remove(score_file)  # Clears file before next run
    os.remove(os.getcwd() + "/ALPHA_flags")  # Remove alpha flag file
    os.remove(os.getcwd() + "/BETA_flags")  # Remove beta flag file
    write_results(results)


# Writing to output file
def write_results(results):
    global verbose
    with open("AB_usage_flex.csv", "w") as f:
        f.write("PDB,ALPHA,BETA\n")
        if verbose:
            print("PDB\tALPHA\tBETA")
        for pdb_id in results:
            f.write(pdb_id + "," + str(results[pdb_id]["ALPHA"]) + "," + str(results[pdb_id]["BETA"]) + "\n")
            if verbose:  # Writes to stdout if verbose
                print(pdb_id + "\t" + str(results[pdb_id]["ALPHA"]) + "\t" + str(results[pdb_id]["BETA"]))


# Creates flag files for InterfaceAnalyzer
def write_flag(chains):
    keys = ["ALPHA", "BETA"]  # Allows for loop of alpha beta
    for chain in keys:
        chain_id = chains[chain] + "_" + chains["MHC"] + chains["peptide"]  # generates chain_id: ex. AC_D
        with open(chain + "_flags", "w") as f:
            f.write("#specific options for InterfaceAnalyzer\n")
            # f.write("-interface " + chain_id + "\n")  # For specific MHC + Peptide vs a or b
            f.write("-fixedchains " + chains["MHC"] + " " + chains["peptide"] + "\n")  # Holds peptide to MHC
            f.write("-compute_packstat true #Rosetta's packstat calculation (slow)\n")
            f.write("-tracer_data_print false #make a score file instead of stdout\n")
            f.write("-out:file:score_only tcr_" + chain + ".sc  #output file\n")
            f.write("-pack_input false #will not relax the input interface residues\n")  # May want to change
            f.write("-pack_separated false #will also pack monomers to calculated dG bin\n")
            f.write("-add_regular_scores_to_scorefile true #will run the rest of rosetta's score func, score12\n\n")
            f.write("#helpful tweeks\n")
            f.write("-atomic_burial_cutoff 0.01 #This is set to help rosetta identify buried polar atoms properly\n")
            f.write("-sasa_calculator_probe_radius 1.4 #This is the default water probe radius for SASA calculations"\
                    ", sometimes lowering the radius helps rosetta more accurately find buried polar atoms\n")
            f.write("-pose_metrics::interface_cutoff 8.0 #This defines how far away a CBeta atom can be from the"\
                    " other chain to be considered an interface residue")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="Directory containing all PDBs", type=str)
    parser.add_argument("-l", "--linux", help="Changes to linux runnable program", default=False,
                        dest="linux", action="store_true")
    parser.add_argument("-m", "--mac", help="Changes to macOS runnable program", default=True,
                        dest="mac", action="store_true")
    parser.add_argument("-c", "--chains", help="Chains in order of MHC,peptide,alpha,beta", default="ACDE",
                        type=str)
    parser.add_argument("-v", "--verbose", help="Writes output to stdout", default=False, action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    # Runnable on GPU server
    if args.linux:
        program_location = "/mnt/fast/Programs/rosetta_src_2020.08.61146_bundle/main/"\
                "source/bin/InterfaceAnalyzer.mpi.linuxgccrelease"
    # Runnable on Austin's MacBook, change to rosetta dir of local machine
    ##### UPDATE IF NOT AUSTIN #####
    elif args.mac:
        program_location = "/Users/austinseamann/Downloads/rosetta_src_2020.08.61146_bundle/main/source/bin/"\
                "InterfaceAnalyzer.macosclangrelease"
    if args.verbose:
        global verbose
        verbose = True
    # Set TCR chains based on user input
    tcr_chains = {"MHC": args.chains[0].upper(), "peptide": args.chains[1].upper(),
                  "ALPHA": args.chains[2].upper(), "BETA": args.chains[3].upper()}
    run_interface(args.directory, program_location, tcr_chains)


if __name__ == '__main__':
    main()
