# This file is a part of the TRain program
# Author: Austin Seamann & Dario Ghersi
# Version: 0.1
# Last Updated: January 12th, 2022

import argparse
import subprocess
import time
import os
from shutil import copyfile
from PDB_Tools_V3 import PdbTools3


####################
# Global Variables #
####################

version = "linuxgccrelease"
program_dir = os.getcwd()
rosetta_dir = ""
no_mpi = False


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


######################################
#               Rigid                #
######################################
# Method: run_rigid()
# Goal: Run rigid docking with either default or provided parameters from user.
# Input:
#   pdb: location of pdb(s)
#   single_file: boolean value of if a single pdb file or a directory of pdb files
def run_rigid(pdb, run_info, native):
    start_relax = time.time()
    print("Running relax...")
    make_dirs()
    run_relax(pdb, run_info["relax"], run_info["cpu_relax"])
    end_relax = (time.time() - start_relax) / 60.0
    start_dock = time.time()
    print("Relaxing Complete")
    pdb_dock = check_score_relax()
    remove_relax(pdb_dock)
    print("Running dock...")
    run_dock(pdb_dock + ".pdb", run_info["docking"], run_info["cpu_docking"], native)
    end_dock = (time.time() - start_dock) / 60.0
    start_refine = time.time()
    print("Docking Complete")
    pdb_refine = check_score_dock()  # In flexible section
    remove_dock(pdb_refine)  # In flexible section
    print("Running refine...")
    run_refine(pdb_refine + ".pdb", run_info["refine"], run_info["cpu_refine"], native)  # In flexible section
    end_refine = (time.time() - start_refine) / 60.0
    remove_refine(check_score_refine())  # In flexible section
    print("Refine Complete")
    with open("time.txt", "w") as file:
        file.write("Relax: " + str(end_relax) + " min : runs " + str(run_info["relax"]) + "\n")
        file.write("Dock: " + str(end_dock) + " min : runs " + str(run_info["docking"]) + "\n")
        file.write("Refine: " + str(end_refine) + " min : runs " + str(run_info["refine"]) + "\n")
        file.write("Total: " + str(end_relax + end_dock + end_refine) + " min\n")


# Method: make_dirs()
# Goal: Create the subfolders necessary for rigid docking
def make_dirs():
    os.mkdir(program_dir + "/output_files/")
    os.mkdir(program_dir + "/output_files/dock/")
    os.mkdir(program_dir + "/output_files/relax")
    os.mkdir(program_dir + "/output_files/refine")


################
#    Relax     #
################
# Method: run_relax()
# Goal: Runs regular relax
def run_relax(pdb, runs, cpus):
    if not no_mpi:
        dir_relax = rosetta_dir + "/main/source/bin/relax.mpi." + version
    else:
        dir_relax = rosetta_dir + "/main/source/bin/relax." + version
    make_relax_file(pdb, runs, cpus)
    if not no_mpi:
        subprocess.run(["mpirun", "--use-hwthread-cpus", dir_relax, "@flag_input_relax"], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    else:
        subprocess.run([dir_relax, "@flag_input_relax"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


# Method: make_relax_file()
# Goal: Generate relax flag file
def make_relax_file(pdb, runs, cpus):
    with open("flag_input_relax", "w") as relax_file:
        relax_file.write("-in:file:s " + pdb + "\n\n")
        relax_file.write("#SBATCH --ntasks=" + str(cpus) + "\n")
        relax_file.write("-nstruct " + str(runs) + " \n\n")
        relax_file.write("#-relax:constrain_relax_to_start_coords\n")
        relax_file.write("#-relax:ramp_constraints false\n\n")
        relax_file.write("-ex1\n-ex2\n\n")
        relax_file.write("-use_input_sc\n-flip_HNQ\n-no_optH false\n\n")
        relax_file.write("-out:path:all output_files/relax\n")


# Method: remove_relax()
# Goal: Remove redundant relax files to conserve space
def remove_relax(best_pdb):
    for pdb in os.listdir(os.getcwd() + "/output_files/relax/"):
        if pdb.endswith(".pdb"):
            if pdb != best_pdb + ".pdb":
                os.remove(os.getcwd() + "/output_files/relax/" + pdb)


# Method: check_score_relax()
# Goal: Determine best relaxed structure to take to docking step
# TODO: Update to new method.
def check_score_relax():
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
# Method: run_dock()
# Goal: Runs rigid dock
def run_dock(pdb, runs, cpus, native):
    if not no_mpi:
        dir_dock = rosetta_dir + "/main/source/bin/docking_protocol.mpi." + version
    else:
        dir_dock = rosetta_dir + "/main/source/bin/docking_protocol." + version
    make_docking_file(pdb, runs, cpus, native)
    if not no_mpi:
        subprocess.run(["mpirun", "--use-hwthread-cpus", dir_dock, "@flag_local_docking"], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    else:
        subprocess.run([dir_dock, "@flag_local_docking"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


# Method: mack_docking_file()
# Goal: Generate flag file for rigid docking
def make_docking_file(pdb, runs, cpus, native):
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
# Method: run_flexible()
# Goal: Run flexible docking with either default or provided parameters from user.
# Input:
#   pdb: location of pdb(s)
#   single_file: boolean value of if a single pdb file or a directory of pdb files
# TODO: Provide options to not remove files, remove files from relax, add split tcr from multidock
def run_flexible(pdb, run_info, native):
    print("Starting!")
    print("Running relax...")
    start = time.time()
    flex_make_dirs()
    split_tcr = ["python3", "PDB_Tools_V3.py", pdb, "--tcr_split_default"]  # Create tcr.pdb file
    split_pmhc = ["python3", "PDB_Tools_V3.py", pdb, "--pmhc_split"]  # Create pmhc.pdb file
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
    check_rmsd(native, "dock")  # Calculate RMSD CA and all-atom - append to score file
#    remove_dock(check_score_dock())  # TODO: provide option to not remove
    print("Running refine...")
    pdb_refine = check_score_dock()
    run_refine(pdb_refine + ".pdb", run_info["refine"], run_info["cpu_refine"], native)
    best_refine = check_score_refine()
    check_rmsd(native, "refine")
#    remove_refine(best_refine)  # TODO: provide option to not remove
    print(best_refine)
    print("DONE!")
    print(f"Time: {(time.time() - start)/60:.0f} mins")


# Method: flex_make_dirs()
# Goal: Create the subfolders necessary for flexible docking
def flex_make_dirs():
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
# Method: run_pmhc_relax()
# Goal: Relax just the pmhc portion of TCR
def run_pmhc_relax(pmhc, runs, cpus):
    if not no_mpi:
        dir_relax = rosetta_dir + "/main/source/bin/relax.mpi." + version
    else:
        dir_relax = rosetta_dir + "/main/source/bin/relax." + version
    make_pmhc_relax_file(pmhc, runs, cpus)
    if not no_mpi:
        subprocess.run(["mpirun", "--use-hwthread-cpus", dir_relax, "@flag_pmhc_relax"], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    else:
        subprocess.run([dir_relax, "@flag_pmhc_relax"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    check_score_pmhc_relax()


# Method: check_score_pmhc_relax()
# Goal: Checks for best pmhc relaxed structure. Returns name and copies it from output to input.
def check_score_pmhc_relax():
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


# Method: make_pmhc_relax_file()
# Goal: make command file for pmhc relax
def make_pmhc_relax_file(pmhc, runs, cpus):
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
# Method: run_xml_relax()
# Goal: Runs XML relax for tcr
def run_xml_relax(tcr, runs, cpus):
    if not no_mpi:
        dir_script = rosetta_dir + "/main/source/bin/rosetta_scripts.mpi." + version
    else:
        dir_script = rosetta_dir + "/main/source/bin/rosetta_scripts." + version
    make_xml_file()
    make_xml_flag(tcr, runs, cpus)
    if not no_mpi:
        subprocess.run(["mpirun", "--use-hwthread-cpus", dir_script, "@flag_xml_relax"], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    else:
        subprocess.run([dir_script, "@flag_xml_relax"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


# Method: make_xml_file()
# Goal: Creates the xml file for normal relax run of tcrs
def make_xml_file():
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


# Method: make_xml_flag()
# Goal: Creates the flag file for the normal xml relax
def make_xml_flag(tcr, runs, cpus):
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
# Method: run_bb_relax()
# Goal: Runs bb relax for tcr
def run_bb_relax(tcr, runs, cpus):
    if not no_mpi:
        dir_relax = rosetta_dir + "/main/source/bin/relax.mpi." + version
    else:
        dir_relax = rosetta_dir + "/main/source/bin/relax." + version
    make_bb_relax_file(tcr, runs, cpus)
    if not no_mpi:
        subprocess.run(["mpirun", "--use-hwthread-cpus", dir_relax, "@flag_bb_tcr_relax"], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    else:
        subprocess.run([dir_relax, "@flag_bb_tcr_relax"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


# Method: make_bb_relax_file()
# Goal: Creates the flag file for the bb rub
def make_bb_relax_file(tcr, runs, cpus):
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
# Method: run_fast_relax()
# Goal: Runs fast relax for tcr
def run_fast_relax(tcr, runs, cpus):
    dir_relax = rosetta_dir + "/main/source/bin/relax.mpi." + version
    make_fast_relax_file(tcr, runs, cpus)
    subprocess.run(["mpirun", "--use-hwthread-cpus", dir_relax, "@flag_fast_relax"], stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)


# Method: make_fast_relax_file
# Goal: Creates flag file for fast relax of tcr
def make_fast_relax_file(tcr, runs, cpus):
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
# Method: run_prepack()
# Goal: Run prepack after relax protocols
def run_prepack(pdb):
    dir_relax = rosetta_dir + "/main/source/bin/docking_prepack_protocol.mpi." + version
    make_ensemble_files()
    make_prepack_file(pdb)
    subprocess.run(["mpirun", "--use-hwthread-cpus", dir_relax, "@flag_ensemble_prepack"], stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)


# Method: make_ensemble_files()
# Goal: Generates ensemble files from the resulting relaxed filed
# Info: Would add to this if further relaxation is done on pMHCs
def make_ensemble_files():
    global program_dir
    with open("pmhc_ensemblelist", "w") as p:
        for filename in sorted(os.listdir(program_dir + "/input_files/pmhc_ensembles/")):
            if filename.endswith(".pdb"):
                p.write("input_files/pmhc_ensembles/" + filename + "\n")
    with open("tcr_ensemblelist", "w") as t:
        for filename in sorted(os.listdir(program_dir + "/input_files/tcr_ensembles")):
            if filename.endswith(".pdb"):
                t.write("input_files/tcr_ensembles/" + filename + "\n")


# Method: make_prepack_file()
# Goal: Make the flag file for generating the resulting prepacked pdb
def make_prepack_file(pdb):
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
# Method: run_flex_dock()
# Goal: Run flexible docking
def run_flex_dock(pdb, runs, cpus, native):
    dir_dock = rosetta_dir + "/main/source/bin/docking_protocol.mpi." + version
    make_flex_docking_file(pdb, runs, cpus, native)
    subprocess.run(["mpirun", "--use-hwthread-cpus", dir_dock, "@flag_ensemble_docking"], stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)


# Method: make_docking_file()
# Goal: Generate flag file for ensemble docking
def make_flex_docking_file(pdb, runs, cpus, native):
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


# Method: remove_dock()
# Goal: Remove files that are not highest scoring
def remove_dock(best_pdb):
    for pdb in os.listdir(os.getcwd() + "/output_files/dock/"):
        if pdb.endswith(".pdb"):
            if pdb != best_pdb + ".pdb":
                os.remove(os.getcwd() + "/output_files/dock/" + pdb)


# Method: check_score_dock()
# Goal: Returns the best pdb from the score log form docking
# TODO: possibly update with new method
def check_score_dock():
    score_dic = {}
    best_pdb = ""
    score = 100000.0
    for file in os.listdir(os.getcwd() + "/output_files/dock"):
        if file.endswith(".sc"):
            score_file = os.getcwd() + "/output_files/dock/" + file
    with open(score_file, "r") as score_read:
        for line in score_read:
            if not line.__contains__("SEQUENCE") and not line.__contains__("total_score"):
                score_dic[line.split()[-1]] = line.split()
    for pdb in score_dic:
        if float(score_dic[pdb][5]) < score:
            best_pdb = score_dic[pdb][-1]
            score = float(score_dic[pdb][5])
    return best_pdb


################
#    Refine    #
################
# Method: run_refine()
# Goal: Run refinement protocol
def run_refine(pdb, runs, cpus, native):
    dir_dock = rosetta_dir + "/main/source/bin/docking_protocol.mpi." + version
    make_refine_file(pdb, runs, cpus, native)
    subprocess.run(["mpirun", "--use-hwthread-cpus", dir_dock, "@flag_local_refine"], stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)


# Method: make_refine_file()
# Goal: Generate flag file for refinement
def make_refine_file(pdb, runs, cpus, native):
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


# Method: check_score_refine()
# Goal: Check for best scoring file with RE
def check_score_refine():
    for file in os.listdir(os.getcwd() + "/output_files/refine"):
        if file.endswith(".fasc"):
            score_file = "output_files/refine/" + file
    express = "tail -n +3 " + score_file + " | tr -s ' ' | sort -u -k6 -r | head -1"
    temp_best = subprocess.run(express, shell=True, stdout=subprocess.PIPE)
    best_pdb = temp_best.stdout.decode('utf-8').split(' ')[-1][:-1]  # Removing new line character
    return best_pdb


# Method: remove_refine()
# Goal: Remove files that are not highest scoring
def remove_refine(refine_best_pdb):
    for pdb in os.listdir(os.getcwd() + "/output_files/refine/"):
        if pdb.endswith(".pdb"):
            if pdb != refine_best_pdb + ".pdb":
                os.remove(os.getcwd() + "/output_files/refine/" + pdb)


####################
#      Multi       #
####################
# Method: run_multi()
# Goal: Manage runs with submitted directory of pdb files
def run_multi(args):
    main_dir = os.getcwd()
    start = time.time()
    os.mkdir("Runs")
    for pdb in sorted(os.listdir(args.pdb)):
        if pdb.endswith(".pdb"):
            prep_dirs(main_dir, args.pdb, pdb)
            os.chdir(main_dir + "/Runs/" + pdb.split(".")[0] + "/")
            print(pdb)
            # Run Flex auto
            par_run = choose_par(args, ["python3", "TCRcoupler.py", pdb], True)
            print(par_run)
            subprocess.run(par_run)
            os.chdir(main_dir)
            current_time = time.time()
            print(f"{pdb}: {(current_time - start) / 3600:.1f} hrs")
    print(f"Total: {(time.time() - start) / 3600:.1f} hrs")


# Method: choose_par()
# Goal: Provide list of executable parameters for run_multi
def choose_par(args, run_list, multi=False):
    if args.flexible:
        run_list.append("-f")
    else:
        run_list.append("-r")
    if args.mac:
        run_list.append("-m")
    else:
        run_list.append("-l")
    if args.nompi:
        run_list.append("--nompi")
    run_list.extend(["-c", str(args.cores), "-a", str(args.relax), "-d", str(args.docking), "-p", str(args.pmhc),
                     "-x", str(args.xml), "-b", str(args.bb), "-s", str(args.fast), "-e", str(args.refine), "-n",
                     str(args.native)])
    if multi:  # If submitting multiple docks with matching native files
        potential_pdb = run_list[2][:4] + ".pdb"
        if potential_pdb in os.listdir(program_dir + "/" + run_list[-1]):
            run_list[-1] = program_dir + "/" + run_list[-1] + "/" + potential_pdb
        else:
            print("Reference Native Structure Not Found For Structure: " + run_list[2])
    return run_list


# Method: prep_dirs()
# Goal: Prepare individual folders for each pdb run
def prep_dirs(main_dir, folder, pdb):
    new_folder = main_dir + "/Runs/" + pdb.split(".")[0] + "/"
    os.mkdir(new_folder)  # Make pdb run file
    coupler = "TCRcoupler.py"
    pdb_tools = "PDB_Tools_V3.py"
    copyfile(main_dir + "/" + coupler, new_folder + coupler)
    # Copy over pdb file
    copyfile(main_dir + "/" + folder + "/" + pdb, new_folder + pdb)
    # Copy over PDB tools
    copyfile(main_dir + "/" + pdb_tools, new_folder + pdb_tools)
    # Copy over config file
    copyfile(main_dir + "/config.ini", new_folder + "config.ini")


####################
#       RMSD       #
####################
# Method: check_rmsd()
# Goal: Calculate alpha carbon and all-atom RMSD values after superimposing native file to pMHC, then append to .sc
# Input:
#   native: location of native file
#   step: either "dock" or "refine"
# Output:
#   Modified .sc file with additional columns "all_rmsd" and "ca_rmsd"
def check_rmsd(native, step="dock"):
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
            print(tool.rmsd(aligned_native,"ACDE", "ACDE", False, True))
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
        print(len(all_rmsds))
        print(len(ca_rmsds))
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
                print(all_rmsds[position])
                print(ca_rmsds[position])
                print(new_line[-1])
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
    parser.add_argument("-l", "--linux", help="Changes to linux runnable program", default=True,
                        dest='linux', action='store_true')
    parser.add_argument("-m", "--mac", help="Changes to mac runnable program", default=False,
                        dest='mac', action='store_true')
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
    parser.add_argument("-d", "--docking", help="(Both) Number of docking runs performed", default=20000,
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
    return parser.parse_args()


def main():
    args = parse_args()
    # Initializing rosetta folder
    global rosetta_dir
    with open("config.ini", "r") as f1:  # Grab rosetta location
        for line in f1:
            if line[:11] == "rosetta_loc":
                rosetta_dir = line[:-1].split("=")[1][1:-1]
    if os.path.isfile(args.pdb):
        # Initialize version of Linux or Mac release
        global version
        if args.mac:
            version = "macosclangrelease"
        elif args.linux:
            version = "linuxgccrelease"
        global no_mpi
        if args.nompi:
            no_mpi = True
        # Initialize flexible or rigid docking & determine if a batch of pdbs or a single file
        run_info = prep_numbers(args.cores, args.flexible, args.relax, args.docking, args.refine, args.pmhc,
                                args.xml, args.bb, args.fast)
        if args.flexible:
            run_flexible(args.pdb, run_info, args.native)
        elif args.rigid:
            run_rigid(args.pdb, run_info, args.native)
    else:
        run_multi(args)


if __name__ == '__main__':
    main()


