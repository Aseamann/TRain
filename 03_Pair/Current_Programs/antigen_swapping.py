from PDB_Tools_V3 import PdbTools3
import os
import sys
import subprocess


# Runs chimera commands based on submitted cmd file and list of pdbs to open
def run_chimera(cmd_file, open_pdb=[]):
    chimera_dir = "/Applications/Chimera.app/Contents/MacOS/chimera"  # Default location of chimera on mac
    run_list = [chimera_dir, "--nogui", "--silent"]  # Needed base arguments
    for each in open_pdb:  # Fill run_list with list of PDBs to open
        run_list.append(each)
    run_list.append(cmd_file)  # Adds in the cmd_file
    subprocess.run(run_list)  # Runs chimera


# Prepares cmd file and submits to run_chimera for aligning pdb to reference
# and producing a separated protein and peptide pdb.
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


def main():
    tool = PdbTools3()  # Initialize PDB tools
    directory = os.getcwd()
    pdb_folder = directory + '/Trimmed/'  # Default structures
    # Loops through all PDBs in full directory
    for pdb in sorted(os.listdir(pdb_folder)):
        if pdb.endswith(".pdb"):
            tool.set_file_name(pdb_folder + pdb)  # Sets PDB in tool
            print(pdb)
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


if __name__ == "__main__":
    main()
