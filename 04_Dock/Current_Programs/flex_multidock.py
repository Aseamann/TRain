import argparse
import subprocess
import time
import os
from shutil import copyfile


version = ""
pdb_dir = ""
program_dir = os.getcwd()


def run_multi(pdbs):
    global program_dir
    global version
    os.mkdir(program_dir + "/Flex_Runs/")
    start = time.time()
    current_time = time.time()
    for pdb in sorted(os.listdir(pdbs)):
        if pdb.endswith(".pdb"):
            prep_dirs(pdb)  # Turn back on
            os.chdir(program_dir + "/Flex_Runs/" + pdb + "/")
            print(pdb)
            # Set up TCR and pMHC pdb files
            split_tcr = ["python3", "PDB_Tools_V3.py", pdb, "--tcr_split"]
            split_pmhc = ["python3", "PDB_Tools_V3.py", pdb, "--pmhc_split"]
            subprocess.run(split_tcr)
            subprocess.run(split_pmhc)
            # Run Flex auto
            subprocess.run(["python3", "flexauto_rosetta.py", pdb, version])
            os.chdir(program_dir)
            current_time = time.time()
            print(f"{pdb}: {(current_time - start)/3600:.1f} hrs")
    print(f"Total: {(time.time() - start)/3600:.1f} hrs")


def prep_dirs(pdb):
    global pdb_dir
    global program_dir
    os.mkdir(program_dir + "/Flex_Runs/" + pdb + "/")  # Make pdb run file
    # Copy over autodock
    auto = "flexauto_rosetta.py"
    pdb_tools = "PDB_Tools_V3.py"
    copyfile(program_dir + "/" + auto, program_dir + "/Flex_Runs/" + pdb + "/" + auto)
    # Copy over pdb file
    copyfile(pdb_dir + "/" + pdb, program_dir + "/Flex_Runs/" + pdb + "/" + pdb)
    # Copy over PDB tools
    copyfile(program_dir + "/" + pdb_tools, program_dir + "/Flex_Runs/" + pdb + "/" + pdb_tools)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdbs", help="Directory of PDB files for docking", type=str)
    parser.add_argument("-l", "--linux", help="Changes to linux runnable program.", default=True,
            dest="linux", action="store_true")
    parser.add_argument("-m", "--mac", help="Changes to mac runnable program.", default=False,
            dest="mac", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    global version
    if args.mac:
        version = "-m"
    elif args.linux:
        version = "-l"
    global pdb_dir
    pdb_dir = os.getcwd() + "/" + args.pdbs
    run_multi(pdb_dir)


if __name__ == '__main__':
    main()

