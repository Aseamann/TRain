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
rosetta_dir = "/mnt/fast/Programs/rosetta_src_2020.08.61146_bundle"


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


######################################
#               Rigid                #
######################################
# Method: run_rigid()
# Goal: Run rigid docking with either default or provided parameters from user.
# Input:
#   pdb: location of pdb(s)
#   single_file: boolean value of if a single pdb file or a directory of pdb files
def run_rigid(pdb, single_file):
    return "help"


######################################
#              Flexible              #
######################################
# Method: run_flexible()
# Goal: Run flexible docking with either default or provided parameters from user.
# Input:
#   pdb: location of pdb(s)
#   single_file: boolean value of if a single pdb file or a directory of pdb files
# TODO: Provide options to not remove files, remove files from relax, add split tcr from multidock
def run_flexible(pdb, single_file):
    if single_file:
        print("Starting!")
        print("Running relax...")
        start = time.time()
        flex_make_dirs()
        run_pmhc_relax("pmhc.pdb")
        run_xml_relax("tcr.pdb")
        run_bb_relax("tcr.pdb")
        run_fast_relax("tcr.pdb")
        print("Preparing prepack...")
        run_prepack(pdb)
        print("Running docking...")
        run_flex_dock(pdb)
        remove_dock(check_score_dock())  # TODO: provide option to not remove
        print("Running refine...")
        pdb_refine = check_score_dock()
        run_refine(pdb_refine + ".pdb")
        best_refine = check_score_refine()
        remove_refine(best_refine)  # TODO: provide option to not remove
        print(best_refine)
        print("DONE!")
        print(f"Time: {(time.time() - start)/60:.0f} mins")
    else:
        flex_multi_make_dirs(pdb)


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
def run_pmhc_relax(pmhc):
    dir_relax = rosetta_dir + "/main/source/bin/relax.mpi." + version
    make_pmhc_relax_file(pmhc)
    subprocess.run(["mpirun", dir_relax, "@flag_pmhc_relax"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
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
def make_pmhc_relax_file(pmhc):
    with open("flag_pmhc_relax", "w") as relax_file:
        relax_file.write("-in:file:s " + pmhc + "\n\n")
        relax_file.write("#SBATCH --ntasks=" + str(cpu_pmhc_relax) + "\n")
        relax_file.write("-nstruct " + str(relax_pmhc_runs) + " \n\n")
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
def run_xml_relax(tcr):
    dir_script = rosetta_dir + "/main/source/bin/rosetta_scripts.mpi." + version
    make_xml_file()
    make_xml_flag(tcr)
    subprocess.run(["mpirun", dir_script, "@flag_xml_relax"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


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
        f.write("\t\t\trandomselect=\"false\" relaxmode=\"relax\" nsample=\"120\"\n")
        f.write("\t\t\tcartesian_minimize=\"false\" />\n\t</MOVERS>\n")
        f.write("\t<APPLY_TO_POSE>\n\t</APPLY_TO_POSE>\n")
        f.write("\t<PROTOCOLS>\n\t\t<Add mover=\"nma\" />\n\t</PROTOCOLS>\n")
        f.write("\t<OUTPUT scorefxn=\"bn15_cart\" />\n")
        f.write("</ROSETTASCRIPTS>\n")


# Method: make_xml_flag()
# Goal: Creates the flag file for the normal xml relax
def make_xml_flag(tcr):
    with open("flag_xml_relax", "w") as f:
        f.write("-in:file:s " + str(tcr) + "\n\n")
        f.write("#SBATCH --ntasks=" + str(cpu_xml_relax) + "\n")
        f.write("-nstruct " + str(relax_xml_runs) + "\n\n")
        f.write("-parser:protocol nma.xml\n\n")
        f.write("-out:path:all input_files/tcr_ensembles\n")
        f.write("-out:suffix _xml\n")


################
#   bb Relax   #
################
# Method: run_bb_relax()
# Goal: Runs bb relax for tcr
def run_bb_relax(tcr):
    dir_relax = rosetta_dir + "/main/source/bin/relax.mpi." + version
    make_bb_relax_file(tcr)
    subprocess.run(["mpirun", dir_relax, "@flag_bb_tcr_relax"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


# Method: make_bb_relax_file()
# Goal: Creates the flag file for the bb rub
def make_bb_relax_file(tcr):
    with open("flag_bb_tcr_relax", "w") as relax_file:
        relax_file.write("-in:file:s " + str(tcr) + "\n\n")
        relax_file.write("#SBATCH --ntasks=" + str(cpu_bb_relax) + "\n")
        relax_file.write("-nstruct " + str(relax_bb_runs) + " \n\n")
        relax_file.write("-backrub:ntrials 20000\n")
        relax_file.write("-backrub:mc_kt 0.6\n\n")  # TODO: Is this supposed to be mc_kt?
        relax_file.write("-out:path:all input_files/tcr_ensembles/\n")
        relax_file.write("-out:suffix _bb\n")


################
#  fast Relax  #
################
# Method: run_fast_relax()
# Goal: Runs fast relax for tcr
def run_fast_relax(tcr):
    dir_relax = rosetta_dir + "/main/source/bin/relax.mpi." + version
    make_fast_relax_file(tcr)
    subprocess.run(["mpirun", dir_relax, "@flag_fast_relax"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


# Method: make_fast_relax_file
# Goal: Creates flag file for fast relax of tcr
def make_fast_relax_file(tcr):
    with open("flag_fast_relax", "w") as relax_file:
        relax_file.write("-in:file:s " + str(tcr) + "\n\n")
        relax_file.write("#SBATCH --ntasks=" + str(cpu_fast_relax) + "\n")
        relax_file.write("-nstruct " + str(relax_fast_runs) + " \n\n")
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
    subprocess.run(["mpirun", dir_relax, "@flag_ensemble_prepack"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


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
def run_flex_dock(pdb):
    dir_dock = rosetta_dir + "/main/source/bin/docking_protocol.mpi." + version
    make_docking_file(pdb)
    subprocess.run(["mpirun", dir_dock, "@flag_ensemble_docking"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


# Method: make_docking_file()
# Goal: Generate flag file for ensemble docking
def make_docking_file(pdb):
    with open("flag_ensemble_docking", "w") as dock_file:
        dock_file.write("-in:file:s output_files/prepack/" + pdb[:-4] + "_prepack_0001.pdb" + "\n")
        dock_file.write("-unboundrot " + pdb + "\n\n")
        dock_file.write("#SBATCH --ntasks=" + str(cpu_dock) + "\n")
        dock_file.write("-nstruct " + str(docking_runs) + " \n\n")
        dock_file.write("-partners AC_DE\n")
        dock_file.write("-dock_pert 3 8\n\n")
        dock_file.write("-spin\n-detect_disulf true\n-rebuild_disulf true\n\n")
        dock_file.write("-ensemble1 pmhc_ensemblelist\n-ensemble2 tcr_ensemblelist\n\n")
        dock_file.write("-ex1\n-ex2aro\n\n")
        dock_file.write("-out:path:all output_files/dock\n")
        dock_file.write("-out:suffix _ensemble_dock\n")


# Method: remove_dock()
# Goal: Remove files that are not highest scoring
def remove_dock(best_pdb):
    for pdb in os.listdir(os.getcwd() + "/output_files/dock/"):
        if pdb.endswith(".pdb"):
            if pdb != best_pdb + ".pdb":
                os.remove(os.getcwd() + "/output_files/dock/" + pdb)


################
#    Refine    #
################
# Method: run_refine()
# Goal: Run refinement protocol
def run_refine(pdb):
    dir_dock = rosetta_dir + "/main/source/bin/docking_protocol.mpi." + version
    make_refine_file(pdb)
    subprocess.run(["mpirun", dir_dock, "@flag_local_refine"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


# Method: make_refine_file()
# Goal: Generate flag file for refinement
def make_refine_file(pdb):
    with open("flag_local_refine", "w") as refine_file:
        refine_file.write("-in:file:s output_files/dock/" + pdb + "\n")
        refine_file.write("#SBATCH --ntasks=" + str(refine_runs) + "\n")
        refine_file.write("-nstruct " + str(refine_runs) + " \n\n")
        refine_file.write("-docking_local_refine\n")
        refine_file.write("-use_input_sc\n\n")
        refine_file.write("-ex1\n-ex2aro\n\n")
        refine_file.write("-out:file:fullatom\n")
        refine_file.write("-out:path:all output_files/refine\n")
        refine_file.write("-out:suffix _local_refine\n")


# Method: check_score_refine()
# Goal: Check for best scoring file with RE
def check_score_refine():
    score_file = "output_files/refine/score_local_refine.fasc"
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


# TODO: make this capable of running all from a single python file (not counting pdb_tools)
def flex_multi_make_dirs(pdb_folder):
    global program_dir
    os.mkdir(program_dir + "/Flex_bb_Runs/" + pdb + "/")  # Make pdb run file
    # Copy over autodock
    auto = "flexauto__rosetta.py"
    pdb_tools = "PDB_Tools_V3.py"
    copyfile(program_dir + "/" + auto, program_dir + "/Flex_bb_Runs/" + pdb + "/" + auto)
    # Copy over pdb file
    copyfile(pdb_dir + "/" + pdb, program_dir + "/Flex_bb_Runs/" + pdb + "/" + pdb)
    # Copy over PDB tools
    copyfile(program_dir + "/" + pdb_tools, program_dir + "/Flex_bb_Runs/" + pdb + "/" + pdb_tools)


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
            with open(__file__, "w") as f1:
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


