#!/usr/bin/python3

######################################################################
# ModelEngine.py -- A component of TRain                             #
# Copyright: Austin Seamann & Dario Ghersi                           #
# Version: 1.1                                                       #
# Last Updated: November 16th, 2024                                  #
# Goal: Submit TCR alpha and beta chain sequences to local install   #
#       of TCRmodel included in Rosetta Suite of Programs            #
#                                                                    #
# Named arguments: -a --alpha (Location of FASTA file containing     #
#                             alpha chains)                          #
#                  -b --beta (Location of FASTA file containing beta #
#                            chains)                                 #
#                  -r --refine (Model with CDR3 loop refinement)     #
#                  -s --start (Starting PDB in FASTA file (ALPHA))   #
#                  -d --dir (Change standard output director for     #
#                           resulting modeled TCRs)                  #
#                  -v --verbose (Readout output from TCRmodel)       #
#                  -m --multi (Split up load to multiple cores)      #
#                  -c --cpus (Number of CPUs allocated | default=all)#
#                  -u --cutoff (Sequence similarity cutoff for       #
#                               templates | default=100)             #
######################################################################


import importlib.resources as pkg_resources
import argparse
import subprocess
import os
import concurrent.futures
import shutil


#################
#    Methods    #
#################
def read_fastas(alpha_in, beta_in, start_pdb, mhc_in=None, peptide_in=None):
    """
    Read in alpha, beta, mhc, and peptide fasta files and generate tcr_seqs dictionary

    Parameters
    __________
    alpha_in : str
        alpha chain fasta file
    beta_in : str
        beta chain fasta file
    start_pdb : str
        optional start position in file for modeling
    mhc_in : str, optional
        mhc chain fasta file
    peptide_in : str, optional
        peptide chain fasta file

    Returns
    _______
    tcr_seqs : dict
        dictionary containing reference to each chain, each header, and each sequence
    """
    # Dictionary created (e.g. tcr_seqs[CHAIN][HEADER] = seq)
    list_files = {"ALPHA": alpha_in, "BETA": beta_in}
    if mhc_in:
        list_files["MHC"] = mhc_in
    if peptide_in:
        list_files["PEPTIDE"] = peptide_in

    tcr_seqs = {chain: {} for chain in list_files}
    if start_pdb != "...":  # Starting further down in fasta file
        start = False
    else:
        start = True

    for chain in list_files:
        with open(list_files[chain], "r") as f1:
            header = ""
            for line in f1:
                if not start:  # Searching for first header
                    if line[0] == ">" and line[1:-1] == str(start_pdb):
                        start = True
                if start:  # Won't start until first header is found
                    if line[0] == ">":
                        header = line[:-1].split(">")[1]
                    elif header != "":  # if we've started our first header
                        if header in tcr_seqs[chain].keys():  # If string already started
                            tcr_seqs[chain][header] += line[:-1]
                        else:  # If first line of seq
                            tcr_seqs[chain][header] = line[:-1]
    return tcr_seqs


def run_pars(pars, verbose, header):
    """
    Run the subprocess based on submitted parameters

    Parameters
    __________
    pars : list
        list of parameters in the correct order needed to run tcrmodel
    verbose : boolean
        if true, allow subprocess to output to terminal
    header : str
        current pdb being modeled
    """
    print("Modeling: " + header)
    if verbose:
        subprocess.run(pars)
    else:
        subprocess.run(pars, stdout=subprocess.DEVNULL, stderr=True)
    print("Done: " + header)


def run_tcrmodel(tcr_seqs, rosetta_loc, refine, verbose, multi, cpu_count, cutoff):
    """
    Manage submission to local install of TCRmodel

    Parameters
    __________
    tcr_seqs : dict
        dictionary containing reference to each chain, each header, and each sequence
    rosetta_loc : str
        location of local rosetta scripts
    refine : boolean
        if true, run additional refinement on CDR3 loops
    verbose : boolean
        if true, allow subprocess to output to terminal
    multi : boolean
        if true, submit to TCRmodel using python ProcessPoolExecutor to run on multiple cores
    cpu_count : int
        number of cpu's to allocate for producing models
    cutoff : int
        sequence similarity cutoff for model templates
    """
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
            if cutoff != 100:
                pars += ["-template_identity_cutoff", str(cutoff)]
            if not multi:
                if verbose:
                    subprocess.run(pars)
                else:
                    subprocess.run(pars, stdout=subprocess.DEVNULL, stderr=True)
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

def run_tcrmodel2(tcr_seqs, tcrmodel2_loc, alphafold_db, alphafold_sif, alphafold_ub_sif, singularity_flag, verbose):
    """
    Manage submission to local install of TCRmodel2

    Parameters
    __________
    tcr_seqs : dict
        dictionary containing reference to each chain, each header, and each sequence
    tcrmodel2_loc : str
        location of local tcrmodel2 scripts
    alphafold_db : str
        location of alphafold database
    alphafold_sif : str
        location of alphafold singularity container -- optional if using singularity install of tcrmodel2
    verbose : boolean
        if true, allow subprocess to output to terminal
    """
    count = 0
    header_count = {}  # Count_id: Header
    for header in tcr_seqs["ALPHA"]:
        alpha_chain = tcr_seqs["ALPHA"][header]
        beta_chain = tcr_seqs["BETA"][header]
        current_tmp = "tmp" + str(count) + "_"
        header_count[current_tmp[:-1]] = header
        pars = [
            "python", os.path.join(tcrmodel2_loc, "run_tcrmodel2.py"),
            "--job_id", current_tmp[:-1],
            "--output_dir", os.getcwd(),
            "--tcra_seq", alpha_chain,
            "--tcrb_seq", beta_chain,
            "--ori_db", alphafold_db,
            "--tp_db", "/opt/tcrmodel2/data/databases",
            "--relax_structures", "True"
        ]
        if "MHC" in tcr_seqs or "PEPTIDE" in tcr_seqs:
            # Update run_tcrmodel2.py to run with run_tcrmodel2_ub_tcr.py
            pars[1] = os.path.join(tcrmodel2_loc, "run_tcrmodel2_ub_tcr.py")
        if "MHC" in tcr_seqs:
            mhc_chain = tcr_seqs["MHC"][header]
            pars += ["--mhca_seq", mhc_chain]
        if "PEPTIDE" in tcr_seqs:
            peptide_chain = tcr_seqs["PEPTIDE"][header]
            pars += ["--pep_seq", peptide_chain]
        if singularity_flag:
            pars = [
                "singularity", "run", "--nv", "-B", alphafold_db, alphafold_ub_sif,
                "--job_id", current_tmp[:-1],
                "--output_dir", os.getcwd(),
                "--tcra_seq", alpha_chain,
                "--tcrb_seq", beta_chain,
                "--ori_db", alphafold_db,
                "--tp_db", "/opt/tcrmodel2/data/databases",
                "--relax_structures", "True"
            ]
            if "MHC" in tcr_seqs or "PEPTIDE" in tcr_seqs:
                pars[5] = alphafold_sif
            if "MHC" in tcr_seqs:
                pars += ["--mhca_seq", mhc_chain]
            if "PEPTIDE" in tcr_seqs:
                pars += ["--pep_seq", peptide_chain]
        if verbose:
            subprocess.run(pars)
        else:
            subprocess.run(pars, stdout=subprocess.DEVNULL, stderr=True)
        count += 1  # Add to tmp count
    for pdb in os.listdir(os.getcwd()):
        pdb_file = os.path.join(os.getcwd(), pdb, "ranked_0.pdb")
        if os.path.exists(pdb_file):
            name = header_count[pdb.split("_")[0]]
            header_count.pop(pdb.split("_")[0])
            os.rename(pdb_file, name + '.pdb')
        # Remove temporary directories
        if os.path.isdir(pdb):
            # Fully remove directory
            shutil.rmtree(pdb)
    # Return failed TCR models
    for key in header_count:
        print("Failed: " + header_count[key])
    


####################
#     Controls     #
####################
def parse_args():
    parser = argparse.ArgumentParser("Submit TCR alpha and beta chain sequences to local install of TCRmodel or TCRmodel2")
    parser.add_argument("--tcrmodel2", help="Use TCRmodel2 instead of TCRmodel", action="store_true", default=False)
    parser.add_argument("-sing", "--singularity", help="Use Singularity container for TCRmodel2", action="store_true", default=False)
    parser.add_argument("-a", "--alpha", help="Location of Fasta file containing alpha chains", type=str)
    parser.add_argument("-b", "--beta", help="Location of Fasta file containing beta chains", type=str)
    parser.add_argument("-m", "--mhc", help="[tcrmodel2] Location of Fasta file containing mhc chain", type=str)
    parser.add_argument("-p", "--peptide", help="[tcrmodel2] Location of Fasta file containing peptide chain", type=str)
    parser.add_argument("-r", "--refine", help="Model with CDR3 loop refinement", action="store_true", default=False)
    parser.add_argument("-s", "--start", help="Starting PDB in Fasta file (Alpha)", type=str, default="...")
    parser.add_argument("-d", "--dir", help="Change standard output directory for resulting modeled TCRs",
                        default="Modeled")
    parser.add_argument("-v", "--verbose", help="Readout output from TCRmodel", default=False, action="store_true")
    parser.add_argument("-c", "--cpus", help="Number of CPUs allocated | default=all", type=int, default=1)
    parser.add_argument("-u", "--cutoff", help="Sequence similarity cutoff for templates | default=100", type=int,
                        default=100)
    return parser.parse_args()


def main():
    args = parse_args()
    # Grab Rosetta and TCRmodel2 install location
    rosetta_dir = ""
    with pkg_resources.path('data', 'config.ini') as p:
        config_file = p
    with open(config_file, "r") as f1:
        for line in f1:
            if "rosetta_loc" in line:
                rosetta_dir = line[:-1].split("=")[1][1:-1]
            if args.tcrmodel2:
                if "tcrmodel2" in line:
                    tcrmodel2_dir = line[:-1].split("=")[1][1:-1]
                if "alphafold_db" in line:
                    alphafold_db = line[:-1].split("=")[1][1:-1]
                if "alphafold_sif" in line:
                    alphafold_sif = line[:-1].split("=")[1][1:-1]
                if "alphafold_ub_sif" in line:
                    alphafold_ub_sif = line[:-1].split("=")[1][1:-1]
    if not rosetta_dir.endswith("/"):
        rosetta_dir += "/"
    # Logic to determine which rosetta binaries are available - prefer MPI version
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
    # Read in FASTA files
    if not args.tcrmodel2:
        # TCRmodel
        tcr_seqs = read_fastas(args.alpha, args.beta, args.start)
        new_dir = os.getcwd() + "/" + args.dir
        if not os.path.exists(new_dir):
            os.mkdir(new_dir)
        os.chdir(new_dir)
        # If CPU > 1, submit to multiple cores
        multi = False
        if args.cpus > 1:
            multi = True
        run_tcrmodel(tcr_seqs, rosetta_dir, args.refine, args.verbose, multi, args.cpus, args.cutoff)
    else:
        # TCRmodel2
        if args.mhc:
            tcr_seqs = read_fastas(args.alpha, args.beta, args.start, args.mhc, args.peptide)
        else:
            tcr_seqs = read_fastas(args.alpha, args.beta, args.start)
        new_dir = os.getcwd() + "/" + args.dir
        if not os.path.exists(new_dir):
            os.mkdir(new_dir)
        os.chdir(new_dir)
        run_tcrmodel2(tcr_seqs, tcrmodel2_dir, alphafold_db, alphafold_sif, alphafold_ub_sif, args.singularity, args.verbose)


if __name__ == '__main__':
    main()
