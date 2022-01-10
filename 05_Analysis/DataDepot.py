# This file is a part of the TRain program
# Author: Austin Seamann & Dario Ghersi
# Version: 0.1
# Last Updated: November 10th, 2021
import argparse
import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import subprocess
from PDB_Tools_V3 import PdbTools3


#################
#    Global     #
#################
rosetta_dir = ""


#################
#    Methods    #
#################
def box_plot_isc(df_0):
    score_df = pd.read_csv(df_0, index_col=0, sep="\t").sort_values('TCR')
    m1_df = score_df.query("TCR_Ant  == 'M1'")
    print(m1_df)
    m1s = {}
    for tcr in m1_df["TCR"].unique():  # for each with M1 TC
        m1s[tcr] = m1_df.query("TCR == '" + tcr + "'")
    for each in m1s:
        current = m1s[each].query("pMHC == '" + each + "'")  # pulls current TCR from table
        m1s[each] = m1s[each].query("pMHC != '" + each + "'")
        ax = sns.boxplot(y=m1s[each]["I_sc"], x=m1s[each]["TCR"])  # Produces box plot
        plt.plot([current["I_sc"], current["I_sc"]], linewidth=2, color="red")
        plt.show()


#################
#   AB usage    #
#################
# Global
verbose = False


# Manages program
def run_interface(tcr_dir, program, chains):
    global verbose
    tool = PdbTools3()  # Initialize PDB tools
    write_flag(chains)  # Creates Alpha and Beta flag files
    results = {}  # Results dictionary that will be used to create output ex. PDBid: {alpha: [scores], beta: [scores]}
    # Score = dG_separated
    # Loops through each PDB in submitted directory
    if os.path.isdir(tcr_dir):
        pdbs = os.listdir(tcr_dir)
    else:
        pdbs = [tcr_dir]
    for pdb in pdbs:
        if pdb.endswith(".pdb"):
            if len(pdbs) > 1:  # If not a single pdb
                file_name = tcr_dir + "/" + pdb
                tool.set_file_name(file_name)
            else:
                file_name = pdb
                tool.set_file_name(file_name)
            results[pdb] = {"ALPHA": -1.0, "BETA": -1.0}
            tool.mute_aa(0, 1000, "B")  # muting the beta chain
            # Alpha
            if verbose:
                if os.path.exists(program):
                    subprocess.run([program, "-s", file_name, "@ALPHA_flags"])  # Runs InterfaceAnalyzer
                else:  # If mpi
                    print("Running with MPI")
                    program_split = program.split(".")[:-1]
                    program_split += ["mpi", program.split(".")[-1]]
                    program = ".".join(program_split)
                    subprocess.run([program, "-s", file_name, "@ALPHA_flags"])  # Runs InterfaceAnalyzer
            else:
                if os.path.exists(program):
                    subprocess.run([program, "-s", file_name, "@ALPHA_flags"],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # Runs InterfaceAnalyzer
                else:  # If mpi
                    print("Running with MPI")
                    program_split = program.split(".")[:-1]
                    program_split += ["mpi", program.split(".")[-1]]
                    program = ".".join(program_split)
                    subprocess.run([program, "-s", file_name, "@ALPHA_flags"],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # Runs InterfaceAnalyzer
            print("Pass A")
            tool.unmute_aa(0, 1000, "B")  # un-muting beta chain
            tool.mute_aa(0, 1000, "A")  # muting the alpha chain
            # Beta
            if verbose:
                if os.path.exists(program):  # If not MPI
                    subprocess.run([program, "-s", file_name, "@BETA_flags"])  # Runs InterfaceAnalyzer
                else:  # If mpi
                    print("Running with MPI")
                    program_split = program.split(".")[:-1]
                    program_split += ["mpi", program.split(".")[-1]]
                    program = ".".join(program_split)
                    subprocess.run([program, "-s", file_name, "@BETA_flags"])  # Runs InterfaceAnalyzer
            else:
                if os.path.exists(program):  # If not MPI
                    subprocess.run([program, "-s", file_name, "@BETA_flags"],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # Runs InterfaceAnalyzer
                else:  # If mpi
                    print("Running with MPI")
                    program_split = program.split(".")[:-1]
                    program_split += ["mpi", program.split(".")[-1]]
                    program = ".".join(program_split)
                    subprocess.run([program, "-s", file_name, "@BETA_flags"],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # Runs InterfaceAnalyzer
            print("Pass B")
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
    with open("AB_Usage.csv", "w") as f:
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
            f.write("-interface " + chain_id + "\n")  # For specific MHC + Peptide vs a or b
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


#################
#     Native    #
#################
def single_graph(score_file, x, y):
    content = subprocess.run("tr -s ' ' < " + score_file + " | tr ' ' ','", shell=True, stdout=subprocess.PIPE)
    if "/" in score_file:  # Detect if not in current dir
        dir_file = "/".join(score_file.split("/")[:-1]) + "/"
    else:
        dir_file = ""
    new_name = score_file.split("/")[-1].split(".")[0] + ".csv"
    with open(dir_file + new_name, "w") as f:
        for line in content.stdout.decode('utf-8').split('\n')[1:]:
            f.write(line + '\n')
    score_df = pd.read_csv(dir_file + new_name, index_col=0)
    sns.relplot(data=score_df, x=x, y=y)
    plt.show()
    os.remove(dir_file + new_name)


def multi_graph(dir_score_files, x, y):
    all_file = os.getcwd() + "/" + dir_score_files + "/all.sc"
    print(all_file)
    header = False
    with open(all_file, "w") as w1:  # Full .sc file
        file_location = os.getcwd() + "/" + dir_score_files  # Location of folder with full path
        for score_file in sorted(os.listdir(file_location)):  # Loop through each score file
            if score_file.endswith(".sc") and score_file != "all.sc":  # Avoid writing over file
                with open(file_location + "/" + score_file, "r") as r1:
                    print(score_file)
                    r1.readline()  # Skip first line of each file
                    if not header:
                        w1.write(r1.readline()[:-1] + "  file\n")
                        header = True
                    for line in r1.readlines()[1:]:
                        # Adding on file name column to all score file
                        w1.write(line[:-1] + "   " + score_file.split(".")[0] + "\n")
    # Convert to csv
    content = subprocess.run("tr -s ' ' < " + dir_score_files + "/all.sc" + " | tr ' ' ','", shell=True,
                             stdout=subprocess.PIPE)
    new_name = all_file.split("/")[-1].split(".")[0] + ".csv"
    with open(dir_score_files + "/" + new_name, "w") as f:
        for line in content.stdout.decode('utf-8').split('\n'):
            f.write(line + '\n')
    score_df = pd.read_csv(dir_score_files + "/" + new_name, index_col=0)
    graph = sns.FacetGrid(score_df, col="file", hue="CAPRI_rank", col_wrap=4, palette="RdYlGn")
    graph.map(sns.scatterplot, x, y)
    graph.add_legend()
    plt.show()
    os.remove(dir_score_files + "/" + "all.sc")
    os.remove(dir_score_files + "/" + "all.csv")


###################
#     Heatmap     #
###################
# Global
first_aa = {}


# Reading energy breakdown excel file
# File: 0: Score, 1: pose_id, 2: resi1, 3: pdbid1, 4: restype1, 5: resi2, 6: pdbid2, 7: restype2, 8: fa_atr
# 9: fa_rep, 10: fa_sol, 11: fa_sol_rep, 12: fa_intra_sol_xover4, 13: ik_ball_wtd, 14: fa_elec, 15: pro_close
# 16: hbond_sr_bb, 17: hbond_lr_bb, 18: hbond_bb_sc, 19: hbond_sc, 20: dslf_fa13, 21: omega, 22: fa_dun
# 23: p_aa_pp, 24: yhh_planarity, 25: ref, 26: rama_prepro, 27: total, 28: description
# Returns: chains = {'A':[['1A','GLY],['2A','SER'],...]}, interactions = {'1A':['1A', '2A',0.0],...],'2A':[..]}
def read_sheet(sheet_in):
    interactions = {}  # Dir of interactions in PDB, Key:
    chains = {}  # Dir of chains in file, Key: Chain_ID, Value: 0 - red_id, 1 - AA
    df = pd.read_csv(sheet_in)
    for key, value in df.iterrows():
        values = value.array
        if values[7] != "onebody":
            if values[3] in interactions.keys():
                interactions[values[3]].append([values[3], values[6], values[27]])
                if values[6] in interactions.keys():
                    interactions[values[6]].append([values[3], values[6], values[27]])
                else:
                    interactions[values[6]] = [[values[3], values[6], values[27]]]
            else:
                interactions[values[3]] = [[values[3], values[6], values[27]]]
                if values[6] in interactions.keys():
                    interactions[values[6]].append([values[3], values[6], values[27]])
                else:
                    interactions[values[6]] = [[values[3], values[6], values[27]]]
        elif values[7] == "onebody":
            if values[3][-1] in chains.keys():
                chains[values[3][-1]].append([values[3], values[4]])
            else:
                first_aa[values[3][-1]] = int(values[3][:-1])
                chains[values[3][-1]] = [[values[3], values[4]]]
    # print(chains)
    # print(interactions)
    return chains, interactions


# Gathers each interaction and produces a dictionary with {'1A':'GLY',...,'105D':'GLY'}
def get_peptide_inter(chains, interactions, chain_in):
    peptide = []
    aa_inter = {}
    for AA_Peptide in chains[chain_in]:
        peptide.append(AA_Peptide)
    for each in peptide:
        if each[0] in aa_inter.keys():
            aa_inter[each[0]].append(interactions[each[0]])
        else:
            aa_inter[each[0]] = [interactions[each[0]]]
    return aa_inter


# Creates the CSV table
def make_peptide_table(aa_list, chain_list, csv_name, chain_in):
    aa_info = {}  # Key: ChainID Value: AA as 3 letter
    cdr_info = {"D": {"CDR1A": range(28, 33), "CDR2A": range(50, 56), "CDR3A": range(91, 100)},
                "E": {"CDR1B": range(29, 33), "CDR2B": range(51, 59), "CDR3B": range(94, 104)}}
    for chain in chain_list:
        for each in chain_list[chain]:
            aa_info[each[0]] = each[1]
    # Produce table
    with open(csv_name, "w") as t1:
        t1.write("Chain: " + chain_in + ",Interaction energy,Partner,CDR Region\n")
        for AA in aa_list:
            for each in aa_list[AA][0]:
                if each[0] == AA:
                    partner = each[1]
                else:
                    partner = each[0]
                if partner[-1] == "D" or partner[-1] == "E":  # Only allows for TCR chains
                    # Writing to file
                    output = aa_info[AA] + " " + str(int(AA[:-1]) - int(first_aa[AA[-1]]) + 1) + "," + str(each[2]) \
                             + "," + aa_info[partner] + " " + partner[:-1] + "," + partner[-1] + "\n"
                    # Writing CDR information
                    for cdr in cdr_info[partner[-1]]:
                        if int(partner[:-1]) in cdr_info[partner[-1]][cdr]:
                            output = output[:-3] + "," + cdr + "\n"
                    t1.write(output)


# Generates data for heatmap based on chain_in and CDR regions
def heatmap_info(aa_list, chain_list, chain_in):
    aa_info = {}  # Key: ChainID Value: aa as 3 letters
    cdr_info = {"D": {"CDR1A": range(28, 33), "CDR2A": range(50, 56), "CDR3A": range(91, 100)},
                "E": {"CDR1B": range(29, 33), "CDR2B": range(51, 59), "CDR3B": range(94, 104)}}
    header_list = []  # Saves the header constructed: numAA
    for chain in chain_list:  # Key: '234E': 'VAL' for every AA
        for each in chain_list[chain]:
            aa_info[each[0]] = each[1]
    file_name = "temp.csv"
    with open(file_name, "w") as t2:
        for chain in cdr_info:  # First line: CDR id's
            for cdr in cdr_info[chain]:
                for num in cdr_info[chain][cdr]:
                    header_list.append(str(num) + chain)
                    t2.write("," + aa_info[str(num) + chain] + " " + str(num))  # Prints header ex. SER 32, VAL 50...
        t2.write("\n")
        for AA in aa_list:  # Loops though each AA in file and writes ones that partner with D or E
            if AA[-1] == chain_in:  # If the interaction is with specified chain
                t2.write(aa_info[AA] + " " + str(AA[:-1]))  # Writes column 0, peptide info
                for tcr_aa in header_list:
                    flag = False  # Flag for if interaction or not
                    for peptide_inter in aa_list[AA][0]:  # Checks each peptide_inter to tcr_aa
                        if peptide_inter[0] == tcr_aa:
                            t2.write("," + str(peptide_inter[2]))  # Writes energy value
                            flag = True
                            break
                        elif peptide_inter[1] == tcr_aa:
                            t2.write("," + str(peptide_inter[2]))  # Writes energy value
                            flag = True
                            break
                    if not flag:
                        t2.write(",0")  # If not broken out of loop, energy 0
                t2.write("\n")
    return file_name


def heatmap(info, chain_in):
    # Read in csv
    df = pd.read_csv(info)
    # Remove columns with only zero values if MHC
    if chain_in == "A":
        # df = df.loc[:, (df != 0).any(axis=0)]
        df = df[(df.sum(axis=1) != 0)]
    # Read in x and y axis labels
    y_axis_labels = list(df.iloc[:, 0])
    x_axis_labels = list(df.iloc[0:, :])[1:]
    # Drop first column
    df = df.iloc[:, 1:]
    # Annotation labels - only when below 0
    labels = df.iloc[:, :].applymap(lambda v: str(v) if v < 0 else '')
    ax = sns.heatmap(df, vmax=0, linewidths=.2, linecolor="grey", xticklabels=x_axis_labels, yticklabels=y_axis_labels,
                     cmap=sns.cubehelix_palette(start=2, rot=0, reverse=True, dark=0, light=1, as_cmap=True),
                     annot=labels, annot_kws={"fontsize": 8, 'rotation': 90}, fmt='')
    plt.xlabel("TCR")
    y_label = "Peptide"
    if chain_in == "A":
        y_label = "MHC"
    plt.ylabel(y_label)
    # Adding in additional tick marks to account for labeling CDR regions
    # if not remove_0:
    ax2 = ax.twiny()
    ax2.set_xlim([0, ax2.get_xlim()[1]])
    ax2.set_xticks(ax.get_xticks())
    x_tick_labels = []
    # Capture first int
    previous_num = int(x_axis_labels[0].split(" ")[-1].split(".")[0])
    cdr_pos = 0
    cdrs = ["CDR1α", "CDR2α", "CDR3α", "CDR1β", "CDR2β", "CDR3β"]
    cdr_color_dic = {"CDR1α": "darkred", "CDR2α": "firebrick", "CDR3α": "red",
                     "CDR1β": "darkblue", "CDR2β": "blue", "CDR3β": "royalblue"}
    colors = []
    results = []
    for each in x_axis_labels:
        num = int(each.split(" ")[-1].split(".")[0])
        if num > previous_num + 1:
            cdr_pos += 1
        if num < previous_num:
            cdr_pos += 1
        previous_num = num
        colors.append(cdr_color_dic[cdrs[cdr_pos]])
        results.append(cdrs[cdr_pos])
        x_tick_labels.append(cdrs[cdr_pos])
    # Set CDR xtick labels
    ax2.set_xticklabels(x_tick_labels, fontsize=8, rotation=90)
    for xtick, color in zip(ax2.get_xticklabels(), colors):
        xtick.set_color(color)
    ax2.tick_params(top=False)
    ax2.tick_params(top=False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    plt.show()


def interface_heatmap(interface_breakdown, mhc):
    # Changes to MHC chain if requested
    chain_in = "C"
    if mhc:
        chain_in = "A"
    chain_list, inter_list = read_sheet(interface_breakdown)
    aa_inter = get_peptide_inter(chain_list, inter_list, chain_in)
    info = heatmap_info(aa_inter, chain_list, chain_in)
    heatmap(info, chain_in)
    os.remove(info)


# Controller method to generate energy breakdown table
def peptide_table(energry_breakdown, mhc):
    # Changes to MHC chain if requested
    chain_in = "C"
    if mhc:
        chain_in = "A"
    chain_list, inter_list = read_sheet(energry_breakdown)
    aa_inter = get_peptide_inter(chain_list, inter_list, chain_in)
    make_peptide_table(aa_inter, chain_list, "output.csv", chain_in)


# Convert tsv to csv
def tsv_to_csv(tsv_in):
    name_in = tsv_in
    # Rename to same name but replace file extension to .csv
    new_name = "".join(tsv_in.split(".")[:-1]) + ".csv"
    # Convert .out (tsv) to .csv
    cmd = "tr -s ' ' < " + name_in + " | tr ' ' ','"
    new_file = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE)
    os.rename(name_in, new_name)
    with open(new_name, "w") as f1:
        f1.write(new_file.stdout.decode('utf-8')[1:])
    return new_name


# Run Rosetta Residue Energy Breakdown program if PDB submitted for heatmap or energy breakdown table
def run_breakdown(pdb_in, mac):
    # Update program location information
    global rosetta_dir
    version = "linuxgccrelease"
    if mac:
        version = "macosclangrelease"
    program_location = "/main/source/bin/residue_energy_breakdown." + version
    program = rosetta_dir + program_location
    # In and out files
    in_file = "-in:file:s " + pdb_in
    file_name = "energy_breakdown.out"
    out_file = "-out:file:silent " + file_name
    # Run process
    try:  # if not mpi
        subprocess.run([program, in_file, out_file], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except:
        print("Running in mpi")
    else:
        program_split = ".".join(program.split(".")[:-1])
        program_split += ["mpi", program.split(".")[-1]]
        program = ".".join(program_split)
        subprocess.run([program, in_file, out_file], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # Convert tsv to csv
    if "/" in pdb_in:
        file_name = "/".join(pdb_in.split("/")[:-1]) + "/" + file_name
    convert_out = tsv_to_csv(file_name)
    return convert_out


####################
#     Controls     #
####################
def parse_args():
    parser = argparse.ArgumentParser()
    # parser.add_argument("-d", "--data", help="TSV containing: pdb_name,TCR,pMHC,Total_Score,Alpha,Beta", type=str)
    # parser.add_argument("-s", "--strip", help="TSV containing: pdb_name,TCR,pMHC,Total_Score,Alpha,Beta")
    # parser.add_argument("-z", "--strip2", help="TSV containing: pdb_name,TCR,pMHC,Total_Score,Alpha,Beta"\
    #                     " | produces strip plot based on alpha/beta score")
    # parser.add_argument("-b", "--box_tcr", help="Box plot of TCR in comparison of native to all others.",
    #                     default=False, action="store_true")
    parser.add_argument("--ab", help="(AB usage) Submit PDB for ab_usage interface scores")
    parser.add_argument("-v", "--verbose", help="(AB usage) Verbose", action="store_true", default=False)
    parser.add_argument("--sc", help="(Native) Score file produced from docking or refinement (or dir of .sc)")
    parser.add_argument("-x", help="(Native) X axis for native structure comparison", type=str)
    parser.add_argument("-y", help="(Native) Y axis for native structure comparison", type=str)
    parser.add_argument("--heatmap", help="(EB) Energy breakdown csv", type=str)
    parser.add_argument("--table", help="(EB) Energy breakdown csv", type=str)
    parser.add_argument("--mhc", help="(EB) Changes energy breakdown to MHC versus peptide", action="store_true"
                        , default=False)
    parser.add_argument("--mac", help="(EB/AB) Changes rosetta from Linux to MacOS version", action="store_true",
                        default=False)
    return parser.parse_args()


def main():
    args = parse_args()
    # Initializing rosetta folder
    global rosetta_dir
    with open("config.ini", "r") as f1:  # Grab rosetta location
        for line in f1:
            if line[:11] == "rosetta_loc":
                rosetta_dir = line[:-1].split("=")[1][1:-1]
    # AB usage output
    if args.ab:
        global verbose
        if args.verbose:
            global verbose
            verbose = True
        version = "linuxgccrelease"
        if args.mac:
            version = "macosclangrelease"
        program_location = rosetta_dir + "/main/source/bin/InterfaceAnalyzer." + version
        # Set TCR chains based on user input
        tcr_chains = {"MHC": "A", "peptide": "C",
                      "ALPHA": "D", "BETA": "E"}
        run_interface(args.ab, program_location, tcr_chains)
    # Native structure comparison
    if args.sc:
        if not os.path.isdir(args.sc):
            single_graph(args.sc, args.x, args.y)
        else:
            multi_graph(args.sc, args.x, args.y)
    # Heatmap or table routing
    if args.heatmap or args.table:
        if args.heatmap:  # Generate heatmap
            breakdown_file = args.heatmap
            if breakdown_file.endswith(".pdb"):  # Run energy breakdown if PDB submitted
                breakdown_file = run_breakdown(args.heatmap, args.mac)
            elif breakdown_file.endswith(".out") or breakdown_file.endswith(".tsv"):  # Convert to csv if tsv
                breakdown_file = tsv_to_csv(breakdown_file)
            interface_heatmap(breakdown_file, args.mhc)
        if args.table:  # Generate breakdown table
            breakdown_file = args.table
            if breakdown_file.endswith(".pdb"):  # Run energy breakdown if PDB submitted
                breakdown_file = run_breakdown(args.table, args.mac)
            elif breakdown_file.endswith(".out") or breakdown_file.endswith(".tsv"):  # Convert to csv if tsv
                breakdown_file = tsv_to_csv(breakdown_file)
            peptide_table(breakdown_file, args.mhc)


if __name__ == '__main__':
    main()
