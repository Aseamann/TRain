#!/usr/bin/python3

######################################################################
# TCRcoupler.py -- A component of TRain                              #
# Copyright: Austin Seamann & Dario Ghersi                           #
# Version: 0.1                                                       #
# Last Updated: February 8th, 2023                                   #
# Goal: Automation for Rosetta Scripts for RosettaDock 3.2 & 4.0     #
#       specific for TCR structures with proper labeling.            #
#                                                                    #
# Positional argument: pdb (PDB file or folder of PDBS for docking)  #
# Named arguments:                                                   #
#                  -q --nompi (Run Rosetta without mpi)              #
#                  -r --rigid (Initializes rigid docking)            #
#                  -f --flexible (Initializes flexible docking)      #
#                  -c --cores (Max number of cpu cores provided)     #
#                  -a --relax ((Rigid) Number of relax runs          #
#                              performed)                            #
#                  -d --docking ((Both) Number of docking runs       #
#                                performed)                          #
#                  -p --pmhc ((Flex) Number of pmhc relax runs)      #
#                  -x --xml ((Flex) Number of xml relax runs or TCR) #
#                  -b --bb ((Flex) Number of backbone rubs runs for  #
#                           TCR)                                     #
#                  -s --fast ((Flex) Number of fast relax runs for   #
#                             TCR)                                   #
#                  -e --refine ((Both) Number of refinement runs     #
#                               done, post dock)                     #
#                  -n --native (Native structure file or folder to   #
#                               run optional comparison to crystal   #
#                               structure; if folder, first 4        #
#                               characters must match)               #
#                  -t --rotation_check                               #
#                     (Check docking permutations to ensure alpha    #
#                      chain is over the N-terminus vs. C-terminus)  #
#                  -u --num_clusters                                 #
#                     (Desired number of clusters (Default: 0, which #
#                      uses an automatic estimate with the eigengap  #
#                      approach)                                     #
#                  -w --cluster (Pairwise spectral cluster top 200   #
#                                I_sc's after docking permutations,  #
#                                suggested for runs with >1000       #
#                                docking runs)                       #
#                  -R --rerun_refine                                 #
#                     (Rerun for refine, helps if needing to adjust  #
#                      number of clusters)                           #
#                  -C --clear_pdbs                                   #
#                     (Remove PDB files to reduce storage demand     #
#                      after runs. Prevents being able to rerun      #
#                      refine with clustering selection of PDBs and  #
#                      rotation check)                               #
######################################################################

import importlib.resources as pkg_resources
import argparse
import subprocess
import time
import os
from shutil import copyfile
import shutil
from util.PDB_Tools_V3 import PdbTools3
from math import sqrt
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy
from sklearn.cluster import SpectralClustering
from sklearn.manifold import MDS
import pickle
import statistics
# For t-SNE clustering
from sklearn.manifold import TSNE
plt.style.use('seaborn-v0_8-whitegrid')


####################
# Global Variables #
####################
program_dir = os.getcwd()
rosetta_dir = ""
no_mpi = False


####################
#     Methods      #
####################
def prep_numbers(cores, flex, relax, docking, refine, pmhc, xml, bb, fast, rotation_check, num_clusters, cluster,
                 rerun_refine, clear_pdbs):
    """
    Determine number of cpus needed for each run, rigid & flexible.

    Parameters
    ----------
    cores : int
        maximum number of cores allotted (both)
    flex : int
        number of flexible runs (rigid)
    relax : int
        number of relax runs (rigid)
    docking : int
        number of docking runs (both)
    refine : int
        number of refinement runs (both)
    pmhc : int
        number of pmhc relax runs (flexible)
    xml : int
        number of xml relax runs (flexible)
    bb : int
        number of backbone rub rubs (flexible)
    fast : int
        number of fast relax runs (flexible)
    rotation_check : boolean
        to check docking permutations to ensure alpha chain is over the N-term vs. C-term
    num_clusters : int
        number of clusters for spectral clustering after tSNE (both)
    cluster : boolean
        if utilizing clustering of docking run results
    rerun_refine : boolean
        only run refinement on selected docked structures
    clear_pdbs : boolean
        don't save pdbs after runs unless best structure

    Returns
    -------
    run_info : dict
        dictionary of parameters needed for running docking protocols
    """
    run_info = {"docking": docking}
    # Determines docking core count
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
    # Save docking scoring options
    run_info["rotation_check"] = rotation_check
    run_info["num_clusters"] = num_clusters
    run_info["cluster_par"] = cluster
    run_info["rerun_refine"] = rerun_refine
    run_info["clear_pdbs"] = clear_pdbs
    # Flexible docking numbers
    if flex:
        run_info["pmhc"] = pmhc  # Number of pmhc relax runs
        # Determine number of cpu cores for pmhc relax
        if cores > pmhc:
            run_info["cpu_pmhc"] = pmhc + 1  # Plus 1 for controller core
        else:
            run_info["cpu_pmhc"] = cores  # Use all cores
        run_info["xml"] = xml  # number of xml relax runs
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


def rosetta_binary(program_in):
    """
    Determine what binaries are built for rosetta based on each program

    Parameters
    ----------
    program_in : str
        Location of Rosetta
    """
    # Update program location information
    global rosetta_dir
    with pkg_resources.path('data', 'config.ini') as p:
        config_file = p
    with open(config_file, "r") as f1:
        for line in f1:
            if line[:11] == "rosetta_loc":
                rosetta_dir = line[:-1].split("=")[1][1:-1]
    if not rosetta_dir.endswith("/"):
        rosetta_dir += "/"
    programs = []
    for program in os.listdir(rosetta_dir + "main/source/bin"):
        if program.startswith(program_in):
            programs.append(program)
    for program in sorted(programs, key=len):
        if program.startswith(program_in + ".mpi"):
            return rosetta_dir + "main/source/bin/" + program  # resulting binary location
        elif program.startswith(program_in + ".default"):
            return rosetta_dir + "main/source/bin/" + program  # resulting binary location
    for program in sorted(programs, key=len):
        if program.startswith(program_in):
            return rosetta_dir + "main/source/bin/" + program


######################################
#               Rigid                #
######################################
def run_rigid(pdb, run_info, native):
    """
    Run rigid docking with either default or provided parameters from user.

    Parameters
    ----------
    pdb : str
        location of pdb(s)
    run_info : dict
        information needed for running RosettaDock
    native : str
        location of native structure file(s) if comparing to native structure
    """
    start_relax = time.time()
    end_relax = 0
    end_dock = 0
    if not run_info["rerun_refine"]:  # If rerunning refine based on new set of parameters
        print("Running relax...")
        make_dirs()
        run_relax(pdb, run_info["relax"], run_info["cpu_relax"])
        end_relax = (time.time() - start_relax) / 60.0
        start_dock = time.time()
        print("Relaxing Complete")
        pdb_dock = check_score_relax()
        if run_info["clear_pdbs"]:
            remove_relax(pdb_dock)
        print("Running dock...")
        run_dock(pdb_dock + ".pdb", run_info["docking"], run_info["cpu_docking"], native)
        end_dock = (time.time() - start_dock) / 60.0
        print("Docking Complete")
    start_refine = time.time()
    pdb_refine = check_score_dock(run_info["cluster_par"], run_info["rotation_check"],
                                  run_info["num_clusters"])
    if run_info["clear_pdbs"]:
        remove_dock(pdb_refine)  # In flexible section
    print("Running refine...")
    run_refine(pdb_refine + ".pdb", run_info["refine"], run_info["cpu_refine"], native)  # In flexible section
    end_refine = (time.time() - start_refine) / 60.0
    # Copy refined structure to output_files
    best_refine = check_score_refine()
    copyfile("output_files/refine/" + best_refine + ".pdb", "output_files/" + best_refine + ".pdb")
    if run_info["clear_pdbs"]:
        remove_refine(best_refine)  # In flexible section
    print("Refine Complete")
    with open("time.txt", "w") as file:
        if end_relax != 0:
            file.write("Relax: " + str(end_relax) + " min : runs " + str(run_info["relax"]) + "\n")
        if end_dock != 0:
            file.write("Dock: " + str(end_dock) + " min : runs " + str(run_info["docking"]) + "\n")
        file.write("Refine: " + str(end_refine) + " min : runs " + str(run_info["refine"]) + "\n")
        file.write("Total: " + str(end_relax + end_dock + end_refine) + " min\n")


def make_dirs():
    """
    Create the subdirectories necessary for rigid docking
    """
    os.mkdir(program_dir + "/output_files/")
    os.mkdir(program_dir + "/output_files/dock/")
    os.mkdir(program_dir + "/output_files/relax")
    os.mkdir(program_dir + "/output_files/refine")


################
#    Relax     #
################
def run_relax(pdb, runs, cpus):
    """
    Manages run of rigid relax option

    Parameters
    ----------
    pdb : str
        location of pdb to begin relaxation on
    runs : int
        number of relax runs to perform
    cpus : int
        number of cpu cores to allocate
    """
    dir_relax = rosetta_binary("relax")
    make_relax_file(pdb, runs, cpus)
    if not no_mpi:
        subprocess.run(["mpirun", dir_relax, "@flag_input_relax"], stdout=subprocess.DEVNULL, stderr=True)
    else:
        subprocess.run([dir_relax, "@flag_input_relax"], stdout=subprocess.DEVNULL, stderr=True)


def make_relax_file(pdb, runs, cpus):
    """
    Generate relax flag file

    Parameters
    ----------
    pdb : str
        location of pdb for relaxation
    runs : int
        number of relax runs to perform
    cpus : int
        number of cpu cores to allocate
    """
    with open("flag_input_relax", "w") as relax_file:
        relax_file.write("-in:file:s " + pdb + "\n\n")
        relax_file.write("#SBATCH --ntasks=" + str(cpus) + "\n")
        relax_file.write("-nstruct " + str(runs) + " \n\n")
        relax_file.write("#-relax:constrain_relax_to_start_coords\n")
        relax_file.write("#-relax:ramp_constraints false\n\n")
        relax_file.write("-ex1\n-ex2\n\n")
        relax_file.write("-use_input_sc\n-flip_HNQ\n-no_optH false\n\n")
        relax_file.write("-out:path:all output_files/relax\n")


def remove_relax(best_pdb):
    """
    Remove redundant relax files to conserve space

    Parameters
    ----------
    best_pdb : str
        remove all pdbs except the best scoring structure
    """
    for pdb in os.listdir(os.getcwd() + "/output_files/relax/"):
        if pdb.endswith(".pdb"):
            if pdb != best_pdb + ".pdb":
                os.remove(os.getcwd() + "/output_files/relax/" + pdb)


# TODO: Update to new method.
def check_score_relax():
    """
    Determine best relaxed structure to take to docking step

    Returns
    -------
    best_pdb : str
        name of pdb with best total_score after relaxation
    """
    score_dic = {}
    best_pdb = ""
    score = 100000.0
    with open("output_files/relax/score.sc", "r") as score_read:
        for line in score_read:
            if not line.__contains__("SEQUENCE") and not line.__contains__("total_score"):
                score_dic[line.split()[-1]] = line.split()
    for pdb in score_dic:
        if float(score_dic[pdb][1]) < score:
            best_pdb = score_dic[pdb][-1]
            score = float(score_dic[pdb][1])
    return best_pdb


################
#     Dock     #
################
def run_dock(pdb, runs, cpus, native):
    """
    Manages docking runs of rigid docking

    Parameters
    ----------
    pdb : str
        location of relaxed pdb to perform docking on
    runs : int
        number of docking permutations to perform
    cpus : int
        number of cpu cores to allocate
    native : str
        location of native structure to calculate RMSD against
    """
    dir_dock = rosetta_binary("docking_protocol")
    make_docking_file(pdb, runs, cpus, native)
    if not no_mpi:
        subprocess.run(["mpirun", dir_dock, "@flag_local_docking"], stdout=subprocess.DEVNULL, stderr=True)
    else:
        subprocess.run([dir_dock, "@flag_local_docking"], stdout=subprocess.DEVNULL, stderr=True)


def make_docking_file(pdb, runs, cpus, native):
    """
    Generate docking flag file

    Parameters
    ----------
    pdb : str
        location of pdb for docking runs
    runs : int
        number of docking permutations to perform
    cpus : int
        number of cpu cores to allocate
    native : str
        location of native structure to calculate RMSD against

    Returns
    -------

    """
    with open("flag_local_docking", "w") as dock_file:
        dock_file.write("-in:file:s output_files/relax/" + pdb + "\n")
        if native != "...":
            dock_file.write("-in:file:native " + native + "\n")
        dock_file.write("-unboundrot output_files/relax/" + pdb + "\n\n")
        dock_file.write("#SBATCH --ntasks=" + str(cpus) + "\n")
        dock_file.write("-nstruct " + str(runs) + " \n\n")
        dock_file.write("-partners AC_DE\n")  # TODO: Update this once adding in chain parameters
        dock_file.write("-dock_pert 3 8\n\n")
        dock_file.write("-ex1\n-ex2aro\n\n")
        dock_file.write("-out:path:all output_files/dock\n")
        dock_file.write("-out:suffix _local_dock\n")


######################################
#              Flexible              #
######################################
def run_flexible(pdb, run_info, native):
    """
    Run flexible docking with either default or provided parameters from user.

    Parameters
    ----------
    pdb : str
        location of pdb(s)
    run_info : str
        information needed for running RosettaDock
    native : str
        location of native structure file(s) if comparing to native structure
    """
    print("Starting!")
    start = time.time()
    if not run_info["rerun_refine"]:
        print("Running relax...")
        flex_make_dirs()
        split_tcr = ["PDB_Tools_V3", pdb, "--tcr_split_default"]  # Create tcr.pdb file
        split_pmhc = ["PDB_Tools_V3", pdb, "--pmhc_split"]  # Create pmhc.pdb file
        subprocess.run(split_tcr)
        subprocess.run(split_pmhc)
        run_pmhc_relax("pmhc.pdb", run_info["pmhc"], run_info["cpu_pmhc"])
        run_xml_relax("tcr.pdb", run_info["xml"], run_info["cpu_xml"])
        run_bb_relax("tcr.pdb", run_info["bb"], run_info["cpu_bb"])
        run_fast_relax("tcr.pdb", run_info["fast"], run_info["cpu_fast"])
        print("Preparing prepack...")
        run_prepack(pdb)
        print("Running docking...")
        run_flex_dock(pdb, run_info["docking"], run_info["cpu_docking"], native)
        if native != "...":
            check_rmsd(native, "dock")  # Calculate RMSD CA and all-atom - append to score file
    print("Running refine...")
    pdb_refine = check_score_dock(run_info["cluster_par"], run_info["rotation_check"],
                                  run_info["num_clusters"])
    if run_info["clear_pdbs"]:
        remove_dock(pdb_refine)
    run_refine(pdb_refine + ".pdb", run_info["refine"], run_info["cpu_refine"], native)
    best_refine = check_score_refine()
    # Copy refined structure to output_files
    copyfile("output_files/refine/" + best_refine + ".pdb", "output_files/" + best_refine + ".pdb")
    if native != "...":
        check_rmsd(native, "refine")
    if run_info["clear_pdbs"]:
        remove_refine(best_refine)
    print(best_refine)
    print("DONE!")
    print(f"Time: {(time.time() - start)/60:.0f} mins")


def flex_make_dirs():
    """
    Create the subdirectories necessary for flexible docking
    """
    os.mkdir(program_dir + "/input_files/")
    os.mkdir(program_dir + "/output_files/")
    os.mkdir(program_dir + "/output_files/dock/")
    os.mkdir(program_dir + "/output_files/relax/")
    os.mkdir(program_dir + "/output_files/refine/")
    os.mkdir(program_dir + "/output_files/prepack/")
    os.mkdir(program_dir + "/input_files/pmhc_ensembles/")
    os.mkdir(program_dir + "/input_files/tcr_ensembles/")


################
#  pMHC Relax  #
################
def run_pmhc_relax(pmhc, runs, cpus):
    """
    Manges run of pMHC relax

    Parameters
    ----------
    pmhc : str
        location of pMHC PDB
    runs : int
        number of relax runs to perform
    cpus : int
        number of cpu cores to allocate
    """
    dir_relax = rosetta_binary("relax")
    make_pmhc_relax_file(pmhc, runs, cpus)
    if not no_mpi:
        subprocess.run(["mpirun", dir_relax, "@flag_pmhc_relax"], stdout=subprocess.DEVNULL, stderr=True)
    else:
        subprocess.run([dir_relax, "@flag_pmhc_relax"], stdout=subprocess.DEVNULL, stderr=True)
    check_score_pmhc_relax()


def check_score_pmhc_relax():
    """
    Checks for best pmhc relaxed structure. Returns name and copies it from output to input.

    Returns
    -------
    best_pdb : str
        name of the best scoring relaxed pMHC structure
    """
    score_dic = {}
    best_pdb = ""
    score = 100000.0
    with open("output_files/relax/score.sc", "r") as score_read:
        for line in score_read:
            if not line.__contains__("SEQUENCE") and not line.__contains__("total_score"):
                score_dic[line.split()[-1]] = line.split()
    for pdb in score_dic:
        if float(score_dic[pdb][1]) < score:
            best_pdb = score_dic[pdb][-1]
            score = float(score_dic[pdb][1])
    copyfile("output_files/relax/" + best_pdb + ".pdb", "input_files/pmhc_ensembles/" + best_pdb + ".pdb")
    return best_pdb


def make_pmhc_relax_file(pmhc, runs, cpus):
    """
    Generate pMHC relax flag file

    Parameters
    ----------
    pmhc : str
        location of pMHC PDB for relaxation
    runs : int
        number of relax runs to perform
    cpus : int
        number of cpu cores to allocate
    """
    with open("flag_pmhc_relax", "w") as relax_file:
        relax_file.write("-in:file:s " + pmhc + "\n\n")
        relax_file.write("#SBATCH --ntasks=" + str(cpus) + "\n")
        relax_file.write("-nstruct " + str(runs) + " \n\n")
        relax_file.write("#-relax:constrain_relax_to_start_coords\n")
        relax_file.write("#-relax:ramp_constraints false\n\n")
        relax_file.write("-ex1\n-ex2\n\n")
        relax_file.write("-use_input_sc\n-flip_HNQ\n-no_optH false\n\n")
        relax_file.write("-out:path:all output_files/relax\n")


################
#   XML Relax  #
################
def run_xml_relax(tcr, runs, cpus):
    """
    Manages run of XML relax of TCR

    Parameters
    ----------
    tcr : str
        location of TCR pdb
    runs : int
        number of XML relax runs to perform
    cpus : int
        number of cpu cores to allocate
    """
    dir_script = rosetta_binary("rosetta_scripts")
    make_xml_file()
    make_xml_flag(tcr, runs, cpus)
    if not no_mpi:
        subprocess.run(["mpirun", dir_script, "@flag_xml_relax"], stdout=subprocess.DEVNULL, stderr=True)
    else:
        subprocess.run([dir_script, "@flag_xml_relax"], stdout=subprocess.DEVNULL, stderr=True)


def make_xml_file():
    """
    Generate XML protocol file
    """
    with open("nma.xml", "w") as f:
        f.write("<ROSETTASCRIPTS>\n")
        f.write("\t<SCOREFXNS>\n\t\t<ScoreFunction name=\"bn15_cart\" weights=\"beta_nov15_cart\" />\n")
        f.write("\t</SCOREFXNS>\n\t<RESIDUE_SELECTORS>\n\t</RESIDUE_SELECTORS>\n")
        f.write("\t<TASKOPERATIONS>\n\t</TASKOPERATIONS>\n\t<FILTERS>\n\t</FILTERS>\n")
        f.write("\t<MOVERS>\n\t\t<NormalModeRelax name=\"nma\" cartesian=\"true\" centroid=\"false\"\n")
        f.write("\t\t\tscorefxn=\"bn15_cart\" nmodes=\"5\" mix_modes=\"true\" pertscale=\"1.0\"\n")
        f.write("\t\t\trandomselect=\"false\" relaxmode=\"relax\" nsample=\"1\"\n")
        f.write("\t\t\tcartesian_minimize=\"false\" />\n\t</MOVERS>\n")
        f.write("\t<APPLY_TO_POSE>\n\t</APPLY_TO_POSE>\n")
        f.write("\t<PROTOCOLS>\n\t\t<Add mover=\"nma\" />\n\t</PROTOCOLS>\n")
        f.write("\t<OUTPUT scorefxn=\"bn15_cart\" />\n")
        f.write("</ROSETTASCRIPTS>\n")


def make_xml_flag(tcr, runs, cpus):
    """
    Generate XML relax flag file

    Parameters
    ----------
    tcr : str
        location of TCR PDB file
    runs : int
        number of XML relax runs to perform
    cpus : int
        number of cpu cores to allocate
    """
    with open("flag_xml_relax", "w") as f:
        f.write("-in:file:s " + str(tcr) + "\n\n")
        f.write("#SBATCH --ntasks=" + str(cpus) + "\n")
        f.write("-nstruct " + str(runs) + "\n\n")
        f.write("-parser:protocol nma.xml\n\n")
        f.write("-out:path:all input_files/tcr_ensembles\n")
        f.write("-out:suffix _xml\n")


################
#   bb Relax   #
################
def run_bb_relax(tcr, runs, cpus):
    """
    Manages run of backbone rub relax of TCR

    Parameters
    ----------
    tcr : str
        location of TCR pdb
    runs : int
        number of backbone rub relax runs to perform
    cpus : int
        number of cpu cores to allocate
    """
    dir_relax = rosetta_binary("relax")
    make_bb_relax_file(tcr, runs, cpus)
    if not no_mpi:
        subprocess.run(["mpirun", dir_relax, "@flag_bb_tcr_relax"], stdout=subprocess.DEVNULL, stderr=True)
    else:
        subprocess.run([dir_relax, "@flag_bb_tcr_relax"], stdout=subprocess.DEVNULL, stderr=True)


def make_bb_relax_file(tcr, runs, cpus):
    """
    Generate backbone rub relax flag file

    Parameters
    ----------
    tcr : str
        location of TCR pdb
    runs : int
        number of backbone rub relax runs to perform
    cpus : int
        number of cpu cores to allocate
    """
    with open("flag_bb_tcr_relax", "w") as relax_file:
        relax_file.write("-in:file:s " + str(tcr) + "\n\n")
        relax_file.write("#SBATCH --ntasks=" + str(cpus) + "\n")
        relax_file.write("-nstruct " + str(runs) + " \n\n")
        relax_file.write("-backrub:ntrials 20000\n")
        relax_file.write("-backrub:mc_kt 0.6\n\n")
        relax_file.write("-out:path:all input_files/tcr_ensembles/\n")
        relax_file.write("-out:suffix _bb\n")


################
#  fast Relax  #
################
def run_fast_relax(tcr, runs, cpus):
    """
    Manages run of fast relax of TCR

    Parameters
    ----------
    tcr : str
        location of TCR pdb
    runs : int
        number of fast relax runs to perform
    cpus : int
        number of cpu cores to allocate
    """
    dir_relax = rosetta_binary("relax")
    make_fast_relax_file(tcr, runs, cpus)
    subprocess.run(["mpirun", dir_relax, "@flag_fast_relax"], stdout=subprocess.DEVNULL, stderr=True)


def make_fast_relax_file(tcr, runs, cpus):
    """
    Generate fast relax flag file

    Parameters
    ----------
    tcr : str
        location of TCR pdb
    runs : int
        number of fast relax runs to perform
    cpus : int
        number of cpu cores to allocate
    """
    with open("flag_fast_relax", "w") as relax_file:
        relax_file.write("-in:file:s " + str(tcr) + "\n\n")
        relax_file.write("#SBATCH --ntasks=" + str(cpus) + "\n")
        relax_file.write("-nstruct " + str(runs) + " \n\n")
        relax_file.write("-relax:thorough\n\n")
        relax_file.write("-out:path:all input_files/tcr_ensembles/\n")
        relax_file.write("-out:suffix _fast\n")


################
#   PrePack    #
################
def run_prepack(pdb):
    """
    Manages run of prepack of TCR and pMHC relax files

    Parameters
    ----------
    pdb : str
        location of complex pre-dock structure PDB
    """
    dir_relax = rosetta_binary("docking_prepack_protocol")
    make_ensemble_files()
    make_prepack_file(pdb)
    subprocess.run(["mpirun", dir_relax, "@flag_ensemble_prepack"], stdout=subprocess.DEVNULL, stderr=True)


def make_ensemble_files():
    """
    Generate ensemble files for the resulting relaxed files
    """
    global program_dir
    with open("pmhc_ensemblelist", "w") as p:
        for filename in sorted(os.listdir(program_dir + "/input_files/pmhc_ensembles/")):
            if filename.endswith(".pdb"):
                p.write("input_files/pmhc_ensembles/" + filename + "\n")
    with open("tcr_ensemblelist", "w") as t:
        for filename in sorted(os.listdir(program_dir + "/input_files/tcr_ensembles")):
            if filename.endswith(".pdb"):
                t.write("input_files/tcr_ensembles/" + filename + "\n")


def make_prepack_file(pdb):
    """
    Generate prepack flag file

    Parameters
    ----------
    pdb : str
        location of complex pre-dock structure PDB
    """
    with open("flag_ensemble_prepack", "w") as f:
        f.write("-in:file:s " + str(pdb) + "\n")
        f.write("-unboundrot " + str(pdb) + "\n\n")
        f.write("-nstruct 1\n")
        f.write("-partners AC_DE\n\n")
        f.write("-ensemble1 pmhc_ensemblelist\n")
        f.write("-ensemble2 tcr_ensemblelist\n")
        f.write("-ex1\n-ex2aro\n\n")
        f.write("-out:path:all output_files/prepack\n")
        f.write("-out:suffix _prepack\n")


################
#     Dock     #
################
def run_flex_dock(pdb, runs, cpus, native):
    """
    Manages run of flexible docking

    Parameters
    ----------
    pdb : str
        location of pdb to begin docking
    runs : int
        number of docking permutations to perform
    cpus : int
        number of cpu cores to allocate
    native : str
        location of native structure for Rosetta to compare to
    """
    dir_dock = rosetta_binary("docking_protocol")
    make_flex_docking_file(pdb, runs, cpus, native)
    subprocess.run(["mpirun", dir_dock, "@flag_ensemble_docking"], stdout=subprocess.DEVNULL, stderr=True)


def make_flex_docking_file(pdb, runs, cpus, native):
    """
    Generate docking flag file

    Parameters
    ----------
    pdb : str
        location of pdb to begin docking
    runs : int
        number of docking permutations to perform
    cpus : int
        number of cpu cores to allocate
    native : str
        location of native structure for Rosetta to compare to
    """
    with open("flag_ensemble_docking", "w") as dock_file:
        dock_file.write("-in:file:s output_files/prepack/" + pdb[:-4] + "_prepack_0001.pdb" + "\n")
        if native != "...":
            dock_file.write("-in:file:native " + native + "\n")
        dock_file.write("-unboundrot " + pdb + "\n\n")
        dock_file.write("#SBATCH --ntasks=" + str(cpus) + "\n")
        dock_file.write("-nstruct " + str(runs) + " \n\n")
        dock_file.write("-partners AC_DE\n")
        dock_file.write("-dock_pert 3 8\n\n")
        dock_file.write("-spin\n-detect_disulf true\n-rebuild_disulf true\n\n")
        dock_file.write("-ensemble1 pmhc_ensemblelist\n-ensemble2 tcr_ensemblelist\n\n")
        dock_file.write("-ex1\n-ex2aro\n\n")
        dock_file.write("-docking_low_res_score motif_dock_score\n")
        dock_file.write("-mh:path:scores_BB_BB " + rosetta_dir
                        + "/main/database/additional_protocol_data/motif_dock/xh_16_\n")
        dock_file.write("-mh:score:use_ss1 false\n")
        dock_file.write("-mh:score:use_ss2 false\n")
        dock_file.write("-mh:score:use_aa1 true\n")
        dock_file.write("-mh:score:use_aa2 true\n")
        dock_file.write("-out:path:all output_files/dock\n")
        dock_file.write("-out:suffix _ensemble_dock\n")


def remove_dock(best_pdb):
    """
    Remove redundant docking permutations to conserve space

    Parameters
    ----------
    best_pdb : str
        remove all pdbs except that best scoring structure
    """
    for pdb in os.listdir(os.getcwd() + "/output_files/dock/"):
        if pdb.endswith(".pdb"):
            if pdb != best_pdb + ".pdb":
                os.remove(os.getcwd() + "/output_files/dock/" + pdb)


def check_score_dock(cluster_par=False, rotation_check=False, num_clusters=0):
    """
    Determine best docked structure to take to relaxation step. Does so with scoring method highlighted
    in manual.

    Parameters
    ----------
    cluster_par : boolean
        optional parameter to perform or not to perform clustering on top 200 structures
    rotation_check : boolean
        optional parameter to perform or not to perform rotation check on top 200 structures
    num_clusters : int
        number of spectral clusters to form on top 200 structures. default of zero will utilize
        eigen gap heuristic to determine number of clusters

    Returns
    -------
    best_pdb : str
        location of pdb selected by scoring method utilized
    """
    dock_dir = os.getcwd() + "/output_files/dock"
    # Produce ordered list of scores
    ordered = check_score(dock_dir)  # All poses {pdb_id: I_sc} in sorted order

    if not cluster_par:
        pickle_name = "dock_rmsds.pickle"  # Will save to pickle and default uses pickle if detected in folder.

        # Return the list of the 200 best scoring structures with confirmed orientation of TCRpMHC
        passed = list(ordered.keys())[:200]  # if not performing rotation check
        # If performing rotation check
        if pickle_name not in os.listdir() and not rotation_check:
            passed = confirm_rotation(dock_dir, ordered)

        # Calculate pairwise rmsd distance matrix
        # Only perform if rotation_check is not selected
        if pickle_name not in os.listdir():
            rmsd_matrix = pairwise_rmsd(dock_dir, passed)
            with open(pickle_name, 'wb') as f:
                pickle.dump(rmsd_matrix, f)
        else:
            with open(pickle_name, 'rb') as f:
                rmsd_matrix = pickle.load(f)

        # TODO: Remove before publish
        # print("Num of Items in Ordered:" + str(len(ordered.keys())))
        # print("Num of Items in Passed:" + str(len(passed)))
        with open("Ordered.tsv", "w") as f:
            for each in ordered:
                f.write(each + "\t" + ordered[each] + "\n")

        # Cluster with tSNE clustering - with results being clustered by spectral clustering
        # best_pdb = cluster_tsne_only(rmsds, ordered, int(num_clusters))

        # Cluster with MDS to generate report on optional alternative poses. Returns list of top alternative poses
        alternative_pdbs = cluster_mds(rmsd_matrix, passed, ordered, int(num_clusters))
        print(alternative_pdbs)
        # Return best pose
        best_pdb = passed[0]  # Best PDB is already provided from rotation check
    elif not rotation_check:  # If not clustering but still checking for rotation
        # Return the list of the 200 best scoring structures with confirmed orientation of TCRpMHC
        passed = confirm_rotation(dock_dir, ordered)
        best_pdb = passed[0]
    else:  # If not clustering or checking for rotation
        best_pdb = next(iter(ordered.keys()))
    return best_pdb


def check_score(directory):
    """
    Create dictionary of structures sorted by I_sc value

    Parameters
    ----------
    directory : str
        directory that contains .sc file produced by rosetta

    Returns
    -------
    ordered : dict
        dictionary of structures sorted by I_sc value {pdb_id: I_sc}
    """
    for file in os.listdir(directory):
        if file.endswith(".sc"):
            score_file = directory + "/" + file
    express = "tail -n +3 " + score_file + " | tr -s ' ' | sort -g -k6"  # k6 == I_sc
    temp_best = subprocess.run(express, shell=True, stdout=subprocess.PIPE)
    sorted_pdbs = temp_best.stdout.decode('utf-8').split("\n")[:-1]
    ordered = {}  # list of {pdb: score}
    for line in sorted_pdbs:
        pdb = line.split(" ")[-1]
        score = line.split(" ")[5]
        ordered[pdb] = score
    return ordered


def euclidean_distance(atom_1, atom_2):
    """
    Calculate the euclidean distance between two atoms

    Parameters
    ----------
    atom_1 : dict
        dictionary of atom_1's XYZ coords
    atom_2 : dict
        dictionary of atom_2's XYZ coords

    Returns
    -------
    distance : float
        distances between two atoms
    """
    distance = sqrt((atom_2['X'] - atom_1['X']) ** 2 + (atom_2['Y'] - atom_1['Y']) ** 2
                    + (atom_2['Z'] - atom_1['Z']) ** 2)
    return distance


def confirm_rotation(directory, ordered):
    """
    Ensure that the alpha chain is on the side of the N-terminus of the peptide and the beta chain
    is on the side of the C-terminus of the peptide.

    Parameters
    ----------
    directory : str
        directory containing to tcr pdbs to be checked
    ordered : dict
        dictionary of docking permutations sorted by I_sc

    Returns
    -------
    passed_pdbs : list
        first 200 structures that passed rotation check
    """
    tool = PdbTools3()
    passed_pdbs = []
    for pdb in ordered.keys():
        file = directory + "/" + pdb + ".pdb"
        tool.set_file_name(file)
        alpha_atoms = tool.get_atoms_on_chain("D")
        beta_atoms = tool.get_atoms_on_chain("E")
        peptide_atoms = tool.get_atoms_on_chain("C")
        alpha_dis = float('inf')
        beta_dis = float('inf')
        for atom in alpha_atoms:
            dis = euclidean_distance(atom, peptide_atoms[0])
            if dis < alpha_dis:
                alpha_dis = dis
        for atom in beta_atoms:
            dis = euclidean_distance(atom, peptide_atoms[0])
            if dis < beta_dis:
                beta_dis = dis
        if alpha_dis < beta_dis:
            passed_pdbs.append(pdb)
    return passed_pdbs[:200]


def pairwise_rmsd(directory, passed):
    """
    Calculate the pairwise rmsd of all pdbs that passed rotation check. Return as a numpy array

    Parameters
    ----------
    directory : str
        location of tcr pdb files
    passed : list
        first 200 structures that passed rotation check

    Returns
    -------
    pairwise_matrix : np array
        numpy array containing the pairwise rmsd values of all 200 structures
    """
    tool = PdbTools3()

    # Perform calculations - only top half of matrix to avoid redundant calculations
    pairwise_matrix = []  # [[value of pdb1vpdb1, value of pdb1vpdb2, ...], [value of pdb2vpdb2, value of pdb2vpdb3]]
    count = 0
    for pdb_1 in passed:
        pdb_1_file = directory + "/" + pdb_1 + ".pdb"
        tool.set_file_name(pdb_1_file)
        pairwise_matrix.append([])
        if count != 0:  # back populate list with previously calculated values to complete matrix
            for i in range(count):
                pairwise_matrix[count].append(pairwise_matrix[i][count])
        for pdb_2 in passed[count:200]:
            pdb_2_file = directory + "/" + pdb_2 + ".pdb"
            tool.set_file_name(pdb_1_file)
            value = tool.rmsd(pdb_2_file, "DE", "DE", ca=True, mute=True)
            pairwise_matrix[count].append(value)
        count += 1

    # Plot first row - all comparisons to top scoring pose
    plt.clf()
    plt.cla()
    plt.close()
    plt.figure(figsize=(15, 10))
    sns.scatterplot(x=range(0, len(pairwise_matrix[0])), y=pairwise_matrix[0], palette="tab10", s=46)
    plt.ylabel("RMSD from Top Pose")
    plt.xlabel("Top Scoring Poses\n(Top Pose to Pose 200)")
    plt.title("RMSD from Top Scoring Pose")
    plt.savefig("Top_Pose_RMSD.png", format="png")

    # Convert to numpy array and return
    return np.array(pairwise_matrix)


def compute_aff_matrix(rmsds_tsne, nn=8):
    """
    Compute a locally adapted affinity matrix for spectral clustering.
    nn is the parameter for the nearest neighbor

    Parameters
    __________
    rmsds_tsne : array
        numpy array of tsne results of pairwise rmsd matrix
    nn : int
        nearest neighbor value for

    Returns
    _______
    aff_mat : array
        affinity matrix to be used for eigen gap heuristic
    """
    # Compute distance matrix
    dist_mat = scipy.spatial.distance_matrix(rmsds_tsne, rmsds_tsne)

    # Compute the sigmas
    sigmas = []
    for i in range(dist_mat.shape[0]):
        temp = np.copy(dist_mat[i])
        temp.sort()
        sigmas.append(np.median(temp[:nn]))

    # Compute the sigma matrix
    sigma_mat = np.matrix(sigmas)
    # sigma_mat = np.asarray(sigmas)
    sigma_mat = np.matmul(np.transpose(sigma_mat), sigma_mat)

    # Build the kernel function
    aff_mat = np.exp(-np.square(dist_mat) / sigma_mat)

    return aff_mat


def eigen_gap_heuristic(aff_mat, max_k=8):
    """
    Compute the eigengap heuristic to estimate the "optimal" number of
    clusters, following von Luxburg

    Parameters
    __________
    aff_mat : array
        affinity matrix
    max_k : int
        maximum k value for calculating eigengap
    """
    # compute the laplacian
    L = np.identity(np.size(aff_mat, 1)) - aff_mat
    # compute the eigenvalues of the laplacian
    w, v = np.linalg.eig(L)
    # sort the eigenvalues
    w.sort()
    # compute the eigengap
    k = min(max_k, len(w))
    if k == 1:
        return 1
    else:
        gaps = np.diff(w[:k])
        max_g = gaps[0]
        num_cl = 2
        for i in range(1, k - 1):
            if gaps[i] > max_g:
                max_g = gaps[i]
                num_cl = i + 2
    return num_cl


def spectral_cluster(rmsd_tsne, num_clusters):
    """
    Perform spectral clustering of tsne results

    Parameters
    ----------
    rmsd_tsne : array
        numpy array of tsne results of pairwise rmsd matrix
    num_clusters : int
        number of cluster to determine with spectral clustering

    Returns
    -------
    labels : list
        labels correlating to top 200 structures based on spectral clustering results
    """
    # Submit to computer affinity matrix
    affMat = compute_aff_matrix(rmsd_tsne)

    # Compute the number of clusters with the eigengap heuristic
    if num_clusters == 0:
        num_clusters = eigen_gap_heuristic(affMat)

    print("\nNumber of Clusters:", num_clusters)

    # Submit to cluster
    sc = SpectralClustering(n_clusters=int(num_clusters), affinity='precomputed',
                            assign_labels='discretize', random_state=1)
    sc.fit(affMat)
    labels = sc.labels_
    return labels, sc


def cluster_tsne_only(rmsds, ordered, num_clusters):
    """
    Cluster using tSNE clustering and then on reduced points, cluster with spectral clustering

    Parameters
    ----------
    rmsds : dict
        dictionary results from pairwise_rmsd method of top 200 structures
    ordered : dict
        dictionary of structures sorted by I_sc value {pdb_id: I_sc}
    num_clusters : int
        number of clusters to determine provided by the user for the spectral clustering step

    Returns
    -------
    best_pdb : str
        location of best scoring pdb based on the results of clustering and best_cluster_tsne

    Output
    ______
    cluster_results.txt
        Results file providing additional data on clustering results
    """
    # Each PDB in sorted order pulled from rmsds
    pdbs = set([k[0] for k in rmsds.keys()])
    # Convert to rmsd matrix
    rmsd_matrix = []
    row_count = 0
    # Populate rmsd matrix
    for pdb in pdbs:
        rmsd_matrix.append([])
        for x in rmsds:
            if x[0] == pdb or x[1] == pdb:
                rmsd_matrix[row_count].append(rmsds[x])

    # Convert to np array
    rmsd_matrix = np.asarray(rmsd_matrix)
    with open("pdbs.pickle", "wb") as f:
        pickle.dump(pdbs, f)
    with open("rmsd_matrix_cluster.pickle", "wb") as f1:
        pickle.dump(rmsd_matrix, f1)

    # Spectral Clustering
    labels, sc = spectral_cluster(rmsd_matrix, num_clusters)

    # Save ordered to pickle
    with open("ordered.pickle", "rb") as f:
        pickle.dump(ordered, f)

    # Collect additional data on cluster and output to file
    lowest_score = best_cluster(ordered, pdbs, labels, sc)

    # Produce plot of original RMSD pairwise data with embedding and labels of spectral clustering
    plt.clf()
    plt.cla()
    plt.close()
    plt.figure(figsize=(15, 10))
    sns.scatterplot(x=rmsd_matrix[:, 0], y=rmsd_matrix[:, 1], hue=labels, palette="tab10", s=46)
    # plt.scatter(x=rmsd_matrix[0], y=rmsd_matrix[1], color='black', marker="x", s=20)
    plt.title("Labels of Spectral Clustering on First Two Dimensions of Distance Matrix")
    plt.savefig("rmsd_spectral.png", format="png")

    return pdbs[0]


def best_cluster(ordered, uniq_keys, labels, rmsd_tsne):
    """
    Determine the cluster with the lowest I_sc and provide the label and pdb. Provide additional data on each cluster

    Parameters
    ----------
    ordered : dict
        dictionary of structures sorted by I_sc value {pdb_id: I_sc}
    uniq_keys : list
        list of the 200 structures in order
    labels : list
        labels correlating to top 200 structures based on spectral clustering results
    rmsd_tsne : dict
        numpy array of tsne results of pairwise rmsd matrix

    Returns
    -------
    uniq_keys[pos[0]] : str
        best scoring pdb
    rmsd_tsne[pos[0]] : int
        top scoring cluster
    """
    # Create dictionary containing cluster label and corresponding pdbs
    cluster_dict = {label: [] for label in list(set(labels))}  # {cluster num: [pdb, score], [pdb, score] ,...}
    for i in range(len(uniq_keys)):
        cluster_dict[labels[i]].append([uniq_keys[i], ordered[uniq_keys[i]]])
    # Collect scores for each pdb in each cluster
    cluster_score = {label: [] for label in list(set(labels))}
    for cluster in cluster_dict:
        for pdb in cluster_dict[cluster]:
            cluster_score[cluster].append(float(ordered[pdb[0]]))
    print("Cluster Score:", cluster_score)

    # Determine average score for each cluster
    cluster_score_avg = {label: [] for label in list(set(labels))}
    for cluster in cluster_dict:
        cluster_score_avg[cluster] = statistics.mean(cluster_score[cluster])

    # Find the centroid of each cluster
    cluster_coor = {label: [] for label in list(set(labels))}  # {cluster num: [coor 1:, coor 1, [X, Y, Z], ...], ...}
    cluster_center = {label: [] for label in list(set(labels))}  # {cluster num: center}
    for cluster in cluster_dict:
        # Collect the coordinates based on cluster
        for i in range(len(rmsd_tsne)):
            if labels[i] == cluster:
                cluster_coor[cluster].append(list(rmsd_tsne[i]))
        array_cluster_coor = np.array(cluster_coor[cluster])
        center = np.mean(array_cluster_coor, axis=0)
        distances = np.sqrt(((array_cluster_coor - center) ** 2).sum(axis=1))
        centroid = array_cluster_coor[np.argmin(distances)]
        pos = list(np.where((rmsd_tsne == centroid).all(axis=1))[0])
        cluster_center[cluster] = [pos[0], centroid, center]

    with open("cluster_results.txt", "w") as f1:
        f1.write("Lowest I_sc Per Cluster:\n")
        f1.write("\tCluster:\tPDB:\tI_sc:\tCoor:\n")
        lowest_score = {}  # {cluster num: [pdb, score]}
        for cluster in cluster_dict:
            sorted_score = sorted(cluster_dict[cluster], key=lambda x:x[1], reverse=True)
            lowest_score[cluster] = [sorted_score[0]]
            print("Sorted_Score:", sorted_score[0][0])
            print("Sorted_Score 2:", sorted_score[0][1])
            print("RMSD TSNE Result:", list(ordered.keys()).index(sorted_score[0][0]))
            print("Coord:", str(rmsd_tsne[list(ordered.keys()).index(sorted_score[0][0])]))
            f1.write("\t" + str(cluster) + "\t" + sorted_score[0][0] + "\t" + str(sorted_score[0][1])
                     + "\t" + str(rmsd_tsne[list(ordered.keys()).index(sorted_score[0][0])]) + "\n")
        f1.write("Cluster Averages:\n")
        f1.write("\tCluster:\tAvg I_sc:\tCentroid PDB:\tCentroid Coor:\tAbsolute Center:\n")
        for cluster in cluster_score_avg:
            f1.write("\t" + str(cluster) + "\t" + str(cluster_score_avg[cluster])
                     + "\t" + str(list(ordered.keys())[cluster_center[cluster][0]]) +
                     "\t" + str(cluster_center[cluster][1]) + "\t" + str(cluster_center[cluster][2]) + "\n")
        f1.write("Cluster PDBs:\n")
        f1.write("\tCluster:\tPDB:\tI_sc:\n")
        for cluster in cluster_dict:
            for pdb in cluster_dict[cluster]:
                f1.write("\t" + str(cluster) + "\t" + pdb[0] + "\t" + str(pdb[1]) + "\n")
    return lowest_score


def cluster_mds(rmsd_matrix, passed, ordered, num_clusters):
    """
    PUT INFORMATION HERE

    Parameters
    __________
    rmsd_matrix : dict
        dictionary results from pairwise_rmsd method of top x structures
    passed : list
        list of 200 pdb ids in order of score
    ordered : dict
        dictionary of structures sorted by I_sc value {pdb_id: I_sc}
    num_clusters : int
        number of clusters to determine provided by the user for the spectral clustering step

    Returns
    -------
    alternative_pdbs : dict
        best scoring alternative pdb
    """
    alternative_pdbs = {}  # Resulting top alternative poses

    # Perform MDS embedding based on rmsd distance matrix
    embedding = MDS(n_components=2, dissimilarity='precomputed', random_state=52)
    x_transformed = embedding.fit_transform(rmsd_matrix)

    # Perform spectral clustering on resulting reduced data
    labels, sc = spectral_cluster(x_transformed, num_clusters)

    # Plot
    plt.clf()
    plt.cla()
    plt.close()
    plt.figure(figsize=(15,10))
    sns.scatterplot(x=x_transformed[:, 0], y=x_transformed[:, 1], hue=labels, palette="tab10", s=46)
    len_passed = len(rmsd_matrix[0])
    plt.title("Pairwise RMSD Comparison of Top " + str(len_passed) + "\nAfter reduction of rmsd distance matrix with MDS and Clustered with Spectral Clustering")
    plt.savefig("rmsd_mds_spectral.png", format="png")

    # Report top alternative pdbs
    # Pull them into a separate directory
    t_dir = "top"
    if os.path.exists(t_dir):
        shutil.rmtree(t_dir)
    os.mkdir(t_dir)
    for pdb in passed:
        shutil.copy("output_files/dock/" + pdb + ".pdb", "top")
    with open("top/top.csv", "w") as f:
        f.write("PDB\tCluster\n")
        for i in range(len(rmsd_matrix[0])):
            f.write(passed[i] + "\t" + str(labels[i]) + "\t" + str(ordered[passed[i]]) + "\n")
            # Save best I_sc per cluster
            if labels[i] in alternative_pdbs.keys():
                if ordered[passed[i]] < alternative_pdbs[labels[i]][0]:
                    alternative_pdbs[labels[i]] = [passed[i], ordered[passed[i]]]
            else:
                alternative_pdbs[labels[i]] = [passed[i], ordered[passed[i]]]

    return alternative_pdbs



################
#    Refine    #
################
def run_refine(pdb, runs, cpus, native):
    """
    Manages run of refinement of top scoring docking permutation

    Parameters
    ----------
    pdb : str
        location of TCR pdb to perform refinement on
    runs : int
        number of refinement runs to perform
    cpus : int
        number of cpu cores to allocate
    native : str
        optional location of native structure to compare against
    """
    dir_dock = rosetta_binary("docking_protocol")
    make_refine_file(pdb, runs, cpus, native)
    subprocess.run(["mpirun", dir_dock, "@flag_local_refine"], stdout=subprocess.DEVNULL, stderr=True)


# Method: make_refine_file()
# Goal: Generate flag file for refinement
def make_refine_file(pdb, runs, cpus, native):
    """
    Generate refinement flag file

    Parameters
    ----------
    pdb : str
        location of TCR pdb to perform refinement on
    runs : int
        number of refinement steps to perform
    cpus : int
        number of cpu cores to allocate
    native : str
        optional location of native structure to compare against
    """
    with open("flag_local_refine", "w") as refine_file:
        refine_file.write("-in:file:s output_files/dock/" + pdb + "\n")
        if native != "...":
            refine_file.write("-in:file:native " + native + "\n")
        refine_file.write("#SBATCH --ntasks=" + str(runs) + "\n")
        refine_file.write("-nstruct " + str(cpus) + " \n\n")
        refine_file.write("-docking_local_refine\n")
        refine_file.write("-use_input_sc\n\n")
        refine_file.write("-partners AC_DE\n")
        refine_file.write("-ex1\n-ex2aro\n\n")
        refine_file.write("-out:file:fullatom\n")
        refine_file.write("-out:path:all output_files/refine\n")
        refine_file.write("-out:suffix _local_refine\n")


def check_score_refine():
    """
    Determine best refined structure as the final output of RosettaDock

    Returns
    -------
    best_pdb : str
        name of pdb with best I_sc after refinement step
    """
    for file in os.listdir(os.getcwd() + "/output_files/refine"):
        if file.endswith(".fasc"):
            score_file = "output_files/refine/" + file
    express = "tail -n +3 " + score_file + " | tr -s ' ' | sort -u -k6 -r | head -1"
    temp_best = subprocess.run(express, shell=True, stdout=subprocess.PIPE)
    best_pdb = temp_best.stdout.decode('utf-8').split(' ')[-1][:-1]  # Removing new line character
    return best_pdb


def remove_refine(refine_best_pdb):
    """
    Remove redundant refinement files to conserve space

    Parameters
    ----------
    refine_best_pdb : str
        remove all pdbs except the best scoring structure
    """
    for pdb in os.listdir(os.getcwd() + "/output_files/refine/"):
        if pdb.endswith(".pdb"):
            if pdb != refine_best_pdb + ".pdb":
                os.remove(os.getcwd() + "/output_files/refine/" + pdb)


####################
#      Multi       #
####################
def run_multi(args):
    """
    Manages runs with submitted directory of pdb files

    Parameters
    ----------
    args : argparse object
        arguments to pass to each instance of docking
    """
    main_dir = os.getcwd()
    start = time.time()
    os.mkdir("Runs")
    for pdb in sorted(os.listdir(args.pdb)):
        if pdb.endswith(".pdb"):
            prep_dirs(main_dir, args.pdb, pdb)
            os.chdir(main_dir + "/Runs/" + pdb.split(".")[0] + "/")
            print(pdb)
            # Run Flex auto
            par_run = choose_par(args, ["TCRcoupler", pdb], True)
            print(par_run)
            subprocess.run(par_run)
            os.chdir(main_dir)
            current_time = time.time()
            print(f"{pdb}: {(current_time - start) / 3600:.1f} hrs")
    print(f"Total: {(time.time() - start) / 3600:.1f} hrs")


def choose_par(args, run_list, multi=False):
    """
    Provide list of executables parameters for run_multi method

    Parameters
    ----------
    args : argparse object
        arguments to pass to each instance of docking
    run_list : list
        list of arguments to pass to run for each docking run
    multi : boolean
        if submitting each run with a native structure to compare against

    Returns
    -------
    run_list : list
        final list of arguments to pass to run for each docking run
    """
    if args.flexible:
        run_list.append("-f")
    else:
        run_list.append("-r")
    if args.nompi:
        run_list.append("--nompi")
    if args.cluster:
        run_list.append("-w")
    if args.rerun_refine:
        run_list.append("-R")
    if args.clear_pdbs:
        run_list.append("-C")
    if args.rotation_check:
        run_list.append("-t")
    run_list.extend(["-c", str(args.cores), "-a", str(args.relax), "-d", str(args.docking), "-p", str(args.pmhc),
                     "-x", str(args.xml), "-b", str(args.bb), "-s", str(args.fast), "-e", str(args.refine), "-u",
                     str(args.num_clusters), "-n", str(args.native)])
    if multi and args.native != "...":  # If submitting multiple docks with matching native files
        potential_pdb = run_list[2][:4] + ".pdb"
        if potential_pdb in os.listdir(program_dir + "/" + run_list[-1]):
            run_list[-1] = program_dir + "/" + run_list[-1] + "/" + potential_pdb
        else:
            print("Reference Native Structure Not Found For Structure: " + run_list[2])
    return run_list


def prep_dirs(main_dir, folder, pdb):
    """
    Prepare individual folders for each pdb run

    Parameters
    ----------
    main_dir : str
        location of main directory for docking run
    folder : str
        name of folder of each docking run
    pdb : str
        location of pdb to move to folder
    """
    new_folder = main_dir + "/Runs/" + pdb.split(".")[0] + "/"
    os.mkdir(new_folder)  # Make pdb run file
    # Copy over pdb file
    copyfile(main_dir + "/" + folder + "/" + pdb, new_folder + pdb)


####################
#       RMSD       #
####################
def check_rmsd(native, step="dock"):
    """
    Calculate alpha carbon and all-atom RMSD values after superimposing native file to pMHC, then append to .sc

    Parameters
    ----------
    native : str
        location of native file
    step : str
        either "dock" or "refine"
    """
    global program_dir
    score_file = ""
    tool = PdbTools3()  # initialize PDBtools
    ca_rmsds = []  # list of carbon alpha rmsd values in-order
    all_rmsds = []  # list of all-atom rmsd values in-order
    # Loop through each docked structure in order
    location = program_dir + "/output_files/" + step
    for pdb in sorted(os.listdir(location)):
        if pdb.endswith(".pdb"):
            tool.set_file_name(native)  # Set to native file
            tool.superimpose(location + "/" + pdb, "AC", "AC")  # Align mhc and peptide of native file "native + '_aligned.pdb'"
            tool.set_file_name(location + "/" + pdb)  # Set to docked pdb
            aligned_native = native[:-4] + "_aligned.pdb"  # Name/location of aligned native file
            ca_rmsds.append(tool.rmsd(aligned_native, "ACDE", "ACDE", True, True))
            all_rmsds.append(tool.rmsd(aligned_native, "ACDE", "ACDE", False, True))
        if pdb.endswith("sc"):  # Find Score file
            score_file = pdb
    # Read in score_file
    with open(location + "/" + score_file, "r") as f1:
        read_score_file = f1.read()
    # Edit score file
    with open(location + "/" + score_file, "w") as f2:
        count = 0
        position = 0
        for line in read_score_file.split("\n")[:-1]:
            new_line = " ".join(line.split())
            if count == 0:  # first line
                f2.write(line + "\n")
                count += 1
            elif count == 1:  # headers
                new_line = new_line.split(" ")
                # Add new headers
                new_line = new_line[:-1] + ["all_rmsd", "ca_rmsd", new_line[-1]]
                f2.write(" ".join(new_line) + "\n")
                count += 1
            else:
                new_line = new_line.split(" ")
                # Add rmsd values
                new_values = [str(all_rmsds[position]), str(ca_rmsds[position]), new_line[-1]]
                new_line = new_line[:-1] + new_values
                f2.write(" ".join(new_line) + "\n")
                position += 1  # move to next pdb


####################
#     Controls     #
####################
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="PDB file or folder of PDBs for docking", type=str)
    parser.add_argument("-q", "--nompi", help="Run Rosetta without mpi, not recommended for large runs", default=False,
                        action="store_true")
    parser.add_argument("-r", "--rigid", help="Initializes flexible docking", default=True,
                        action="store_true")
    parser.add_argument("-f", "--flexible", help="Initializes flexible docking", default=False,
                        action='store_true')
    parser.add_argument("-c", "--cores", help="Max number of cpu cores provided, default=all",
                        default=os.cpu_count(), type=int)
    parser.add_argument("-a", "--relax", help="(Rigid) Number of relax runs performed", default=100,
                        type=int)
    parser.add_argument("-d", "--docking", help="(Both) Number of docking runs performed", default=10000,
                        type=int)
    parser.add_argument("-p", "--pmhc", help="(Flex) Number of pmhc relax runs", default=100,
                        type=int)
    parser.add_argument("-x", "--xml", help="(Flex) Number of xml relax runs for TCR", default=40,
                        type=int)
    parser.add_argument("-b", "--bb", help="(Flex) Number of back bone rub runs for TCR", default=30,
                        type=int)
    parser.add_argument("-s", "--fast", help="(Flex) Number of fast relax runs for TCR", default=30,
                        type=int)
    parser.add_argument("-e", "--refine", help="(Both) Number of refinement runs done, post dock", default=100,
                        type=int)
    parser.add_argument("-n", "--native", help="Native structure file or folder to run optional comparison to crystal" \
                                               "structure; if folder, first 4 characters must match", type=str,
                                               default="...")
    parser.add_argument("-t", "--rotation_check", help="(Select to not perform) Check docking permutations to " \
                                                       "ensure alpha chain is over the N-terminus vs. the C-terminus",
                        default=False, action='store_true')
    parser.add_argument("-u", "--num_clusters",
                        help="Desired number of clusters (Default: 0, which uses an automatic estimate with the " \
                             "eigengap approach)", default=0)
    parser.add_argument("-w", "--cluster",
                        help="(Select to not perform) Pairwise spectral cluster top x structures I_sc's after " \
                             "docking permutations,suggested for runs with >1000 docking runs", default=False,
                             action='store_true')
    parser.add_argument("-R", "--rerun_refine", help="Rerun from refine, helps if needing to adjust number of clusters",
                        default=False, action='store_true')
    parser.add_argument("-C", "--clear_pdbs", help="Remove PDB files to reduce storage demand after runs. "\
                        "Prevents being able to rerun refine with clustering selection of PDBs and rotation check.",
                        default=False, action='store_true')
    return parser.parse_args()


def main():
    args = parse_args()
    # Initializing rosetta folder
    global rosetta_dir
    with pkg_resources.path('data', 'config.ini') as p:
        config_file = p
    with open(config_file, "r") as f1:  # Grab rosetta location
        for line in f1:
            if line[:11] == "rosetta_loc":
                rosetta_dir = line[:-1].split("=")[1][1:-1]
    print(rosetta_dir)
    if os.path.isfile(args.pdb):
        global no_mpi
        if args.nompi:
            no_mpi = True
        # Initialize flexible or rigid docking & determine if a batch of pdbs or a single file
        run_info = prep_numbers(args.cores, args.flexible, args.relax, args.docking, args.refine, args.pmhc,
                                args.xml, args.bb, args.fast, args.rotation_check, args.num_clusters, args.cluster,
                                args.rerun_refine, args.clear_pdbs)
        if args.flexible:
            run_flexible(args.pdb, run_info, args.native)
        elif args.rigid:
            run_rigid(args.pdb, run_info, args.native)
    else:
        run_multi(args)


if __name__ == '__main__':
    main()



