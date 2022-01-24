# This file is a part of the TRain program
# Author: Austin Seamann & Dario Ghersi
# Version: 0.2
# Last Updated: January 19th, 2022

import argparse
import subprocess
import os
import concurrent.futures


# Read in alpha and beta fasta files and generate tcr_seqs dictionary
def read_fastas(alpha_in, beta_in, start_pdb):
    # Dictionary created (e.g. tcr_seqs[CHAIN][HEADER] = seq)
    list_files = {"ALPHA": alpha_in, "BETA": beta_in}
    tcr_seqs = {"ALPHA": {}, "BETA": {}}
    if start_pdb != "...":  # Starting further down in fasta file
        start = False
    else:
        start = True
    for chain in list_files:
        with open(list_files[chain], "r") as f1:
            header = ""
            for line in f1:
                if not start:  # Searching for first header
                    if line.split(">")[:-1].split(">")[1] == header:
                        start = True
                if start:  # Wont start until first header is found
                    if line[0] == ">":
                        header = line[:-1].split(">")[1]
                    elif header != "":  # if we've started our first header
                        if header in tcr_seqs[chain].keys():  # If string already started
                            tcr_seqs[chain][header] += line[:-1]
                        else:  # If first line of seq
                            tcr_seqs[chain][header] = line[:-1]
    return tcr_seqs


def run_pars(pars, verbose, header):
    print("Modeling: " + header)
    if verbose:
        subprocess.run(pars)
    else:
        subprocess.run(pars, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print("Done: " + header)


def run_tcrmodel(tcr_seqs, rosetta_loc, refine, verbose, multi, cpu_count):
    count = 0
    header_count = {}  # Count_id: Header
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_count) as executor:
        for header in tcr_seqs["ALPHA"]:
            alpha_chain = tcr_seqs["ALPHA"][header]
            beta_chain = tcr_seqs["BETA"][header]
            # Command line parameters
            current_tmp = "tmp" + str(count) + "_"
            header_count[current_tmp[:-1]] = header
            pars = [rosetta_loc, "-alpha", alpha_chain, "-beta", beta_chain, "-out::prefix", current_tmp]
            if refine:  # Additional CDR3 loop remodeling and refinement
                pars += ["-remodel_tcr_cdr3_loops", "-refine_tcr_cdr3_loops"]
            if not multi:
                if verbose:
                    subprocess.run(pars)
                else:
                    subprocess.run(pars, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            else:
                f1 = executor.submit(run_pars, pars, verbose, header)
            count += 1  # Add to tmp count
    for pdb in os.listdir(os.getcwd()):
        if pdb.endswith("tcrmodel.pdb"):
            name = header_count[pdb.split("_")[0]]
            header_count.pop(pdb.split("_")[0])
            os.rename(pdb, name + '.pdb')  # Rename to header name
        elif pdb.startswith("tmp"):
            os.remove(pdb)
    # Return failed TCR models
    for key in header_count:
        print("Failed: " + header_count[key])


####################
#     Controls     #
####################
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--alpha", help="Location of Fasta file containing alpha chains", type=str)
    parser.add_argument("-b", "--beta", help="Location of Fasta file containing beta chains", type=str)
    parser.add_argument("-r", "--refine", help="Model with CDR3 loop refinement", action="store_true", default=False)
    parser.add_argument("-s", "--start", help="Starting PDB in Fasta file (Alpha)", type=str, default="...")
    parser.add_argument("-d", "--dir", help="Change standard output directory for resulting modeled TCRs",
                        default="Modeled")
    parser.add_argument("-v", "--verbose", help="Readout output from TCRmodel", default=False, action="store_true")
    parser.add_argument("-m", "--multi", help="Split up load to multiple cores", action="store_true", default=False)
    parser.add_argument("-c", "--cpus", help="Number of CPUs allocated | default=all", type=int, default=os.cpu_count())
    return parser.parse_args()


def main():
    args = parse_args()
    # Grab Rosetta install location
    rosetta_dir = ""
    with open("config.ini", "r") as f1:
        for line in f1:
            if line[:11] == "rosetta_loc":
                rosetta_dir = line[:-1].split("=")[1][1:-1]
    if not rosetta_dir.endswith("/"):
        rosetta_dir += "/"
    programs = []
    for program in os.listdir(rosetta_dir + "main/source/bin"):
        if program.startswith("tcrmodel."):
            programs.append(program)
    for program in sorted(programs, key=len):
        if program.startswith("tcrmodel.mpi"):
            rosetta_dir += "main/source/bin/" + program
            break
        elif not program.startswith("tcrmodel.default"):
            rosetta_dir += "main/source/bin/" + program
    tcr_seqs = read_fastas(args.alpha, args.beta, args.start)
    new_dir = os.getcwd() + "/" + args.dir
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
    os.chdir(new_dir)
    run_tcrmodel(tcr_seqs, rosetta_dir, args.refine, args.verbose, args.multi, args.cpus)


if __name__ == '__main__':
    main()
