import argparse
import subprocess
import time
import os
from shutil import copyfile


rosetta_dir = ""
program_dir = os.getcwd()
version = ""
relax_runs = 100
docking_runs = 10000
refine_runs = 100

cpu_relax = 101
cpu_dock = 101
cpu_refine = 101


def run_rosetta(pdb):
    start_relax = time.time()
    print("Running relax...")
    make_dirs()
    run_relax(pdb)
    end_relax = (time.time() - start_relax) / 60.0
    start_dock = time.time()
    print("Relaxing Complete")
    pdb_dock = check_score_relax()
    remove_relax(pdb_dock)
    print("Running dock...")
    run_dock(pdb_dock + ".pdb")
    end_dock = (time.time() - start_dock) / 60.0
    start_refine = time.time()
    print("Docking Complete")
    pdb_refine = check_score_dock()
    remove_dock(pdb_refine)
    print("Running refine...")
    run_refine(pdb_refine + ".pdb")
    end_refine = (time.time() - start_refine) / 60.0
    remove_refine(check_score_refine())
    print("Refine Complete")
    with open("time.txt", "w") as file:
        file.write("Relax: " + str(end_relax) + " min : runs " + str(relax_runs) + "\n")
        file.write("Dock: " + str(end_dock) + " min : runs " + str(docking_runs) + "\n")
        file.write("Refine: " + str(end_refine) + " min : runs " + str(refine_runs) + "\n")
        file.write("Total: " + str(end_relax + end_dock + end_refine) + " min\n")

def make_dirs():
    os.mkdir(program_dir + "/output_files/")
    os.mkdir(program_dir + "/output_files/dock/")
    os.mkdir(program_dir + "/output_files/relax")
    os.mkdir(program_dir + "/output_files/refine")

##############################################
# Runs regular relax
def run_relax(pdb):
    dir_relax = rosetta_dir + "/main/source/bin/relax.mpi." + version
    make_relax_file(pdb)
    process = subprocess.run(["mpirun", dir_relax, "@flag_input_relax"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def remove_relax(best_pdb):
    for pdb in os.listdir(os.getcwd() + "/output_files/relax/"):
        if pdb.endswith(".pdb"):
            if pdb != best_pdb + ".pdb":
                os.remove(os.getcwd() + "/output_files/relax/" + pdb)


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


def remove_dock(best_pdb):
    for pdb in os.listdir(os.getcwd() + "/output_files/dock/"):
        if pdb.endswith(".pdb"):
            if pdb != best_pdb + ".pdb":
                os.remove(os.getcwd() + "/output_files/dock/" + pdb)


def check_score_dock():
    score_dic = {}
    best_pdb = ""
    score = 100000.0
    with open("output_files/dock/score_local_dock.sc", "r") as score_read:
        for line in score_read:
            if not line.__contains__("SEQUENCE") and not line.__contains__("total_score"):
                score_dic[line.split()[-1]] = line.split()
    for pdb in score_dic:
        if float(score_dic[pdb][5]) < score:
            best_pdb = score_dic[pdb][-1]
            score = float(score_dic[pdb][5])
    return best_pdb


def run_dock(pdb):
    dir_dock = rosetta_dir + "/main/source/bin/docking_protocol.mpi." + version
    make_docking_file(pdb)
    process = subprocess.run(["mpirun", dir_dock, "@flag_local_docking"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def run_refine(pdb):
    dir_dock = rosetta_dir + "/main/source/bin/docking_protocol.mpi." + version
    make_refine_file(pdb)
    process = subprocess.run(["mpirun", dir_dock, "@flag_local_refine"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def remove_refine(best_pdb):
    for pdb in os.listdir(os.getcwd() + "/output_files/refine/"):
        if pdb.endswith(".pdb"):
            if pdb != best_pdb + ".pdb":
                os.remove(os.getcwd() + "/output_files/refine/" + pdb)


def check_score_refine():
    score_dic = {}
    best_pdb = ""
    score = 100000.0
    with open("output_files/refine/score_local_refine.fasc", "r") as score_read:
        for line in score_read:
            if not line.__contains__("SEQUENCE") and not line.__contains__("total_score"):
                score_dic[line.split()[-1]] = line.split()
    for pdb in score_dic:
        if float(score_dic[pdb][5]) < score:
            best_pdb = score_dic[pdb][-1]
            score = float(score_dic[pdb][5])
    return best_pdb


def make_relax_file(pdb):
    with open("flag_input_relax", "w") as relax_file:
        relax_file.write("-in:file:s " + pdb + "\n\n")
        relax_file.write("#SBATCH --ntasks=" + str(cpu_relax) + "\n")
        relax_file.write("-nstruct " + str(relax_runs) + " \n\n")
        relax_file.write("#-relax:constrain_relax_to_start_coords\n")
        relax_file.write("#-relax:ramp_constraints false\n\n")
        relax_file.write("-ex1\n-ex2\n\n")
        relax_file.write("-use_input_sc\n-flip_HNQ\n-no_optH false\n\n")
        relax_file.write("-out:path:all output_files/relax\n")


def make_docking_file(pdb):
    with open("flag_local_docking", "w") as dock_file:
        dock_file.write("-in:file:s output_files/relax/" + pdb + "\n")
        dock_file.write("-unboundrot output_files/relax/" + pdb + "\n\n")
        dock_file.write("#SBATCH --ntasks=" + str(cpu_dock) + "\n")
        dock_file.write("-nstruct " + str(docking_runs) + " \n\n")
        dock_file.write("-partners DE_AC\n")
        dock_file.write("-dock_pert 3 8\n\n")
        dock_file.write("-ex1\n-ex2aro\n\n")
        dock_file.write("-out:path:all output_files/dock\n")
        dock_file.write("-out:suffix _local_dock\n")


def make_refine_file(pdb):
    with open("flag_local_refine", "w") as refine_file:
        refine_file.write("-in:file:s output_files/dock/" + pdb + "\n")
        refine_file.write("#SBATCH --ntasks=" + str(cpu_refine) + "\n")
        refine_file.write("-nstruct " + str(refine_runs) + " \n\n")
        refine_file.write("-docking_local_refine\n")
        refine_file.write("-use_input_sc\n\n")
        refine_file.write("-ex1\n-ex2aro\n\n")
        refine_file.write("-out:file:fullatom\n")
        refine_file.write("-out:path:all output_files/refine\n")
        refine_file.write("-out:suffix _local_refine\n")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="PDB file of TCRpMHC", type=str)
    parser.add_argument("-l", "--linux", help="Changes to linux runnable program.", default=True,
                        dest='linux', action='store_true')
    parser.add_argument("-m", "--mac", help="Changes to mac runnable program.", default=False,
                        dest="mac", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    global version
    if args.mac:
        version = "macosclangrelease"
    elif args.linux:
        version = "linuxgccrelease"
    global rosetta_dir
    rosetta_dir = "/mnt/fast/Programs/rosetta_src_2020.08.61146_bundle"
    run_rosetta(args.pdb)


if __name__ == '__main__':
    main()

