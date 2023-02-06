#!/usr/bin/python3

######################################################################
# DataDepot.py -- A component of TRain                               #
# Copyright: Austin Seamann & Dario Ghersi                           #
# Version: 0.1                                                       #
# Last Updated: April 3rd, 2022                                      #
# Goal: Perform analysis of TCR PDB files and data from docking      #
#                                                                    #
# Named arguments: -a --ab ((AB usage) Submit PDB for ab_usage       #
#                           interface scores)                        #
#                  -v --verbose ((AB usage) Verbose)                 #
#                  -s --sc ((Native) Score file produced from        #
#                           docking or refinement (or dir of .sc))   #
#                  -x --xaxis ((Native) X axis for native structure  #
#                              comparison)                           #
#                  -y --yaxis ((Native) Y axis for native structure  #
#                              comparison)                           #
#                  --ymax ((Native) Y axis value maximum)            #
#                  --ymin ((Native) Y axis value minimum)            #
#                  --xmax ((Native) X axis value maximum)            #
#                  --xmin ((Native) X axis value minimum)            #
#                  -p --heatmap_dist ((DHM) Distance Breakdown TCR   #
#                                     PDB File)                      #
#                  -d --distance ((DHM) Distance Breakdown TCR PDB   #
#                                 File)                              #
#                  -c --alpha_carbon ((DHM) Alpha carbon only)       #
#                  -e --heatmap_energy ((EB) TCR PDB File)           #
#                  -t --table ((EB) TCR PDB File)                    #
#                  -m --mhc ((DHM/EB) Changes energy breakdown to    #
#                            MHC versus peptide)                     #
#                  -f --fontsize ((DHM/EB) Adjust font size)         #
######################################################################

import argparse
import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import subprocess
from util.PDB_Tools_V3 import PdbTools3


#################
#    Global     #
#################
rosetta_dir = ""
verbose = False


#################
#    Methods    #
#################
#################
#   AB usage    #
#################
def run_interface(tcr_dir, chains):
    """
    Manges the run of interface breakdown

    Parameters
    ----------
    tcr_dir : str
        directory of where the TCR files are
    chains : str
        chains to compare
    """
    global rosetta_dir
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
            tool.mute_aa(0, 100000, "E")  # muting the beta chain
            # Alpha
            if verbose:
                subprocess.run([rosetta_dir, "-s", file_name, "@ALPHA_flags"])  # Runs InterfaceAnalyzer
            else:
                subprocess.run([rosetta_dir, "-s", file_name, "@ALPHA_flags"],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # Runs InterfaceAnalyzer
            tool.unmute_aa(0, 100000, "E")  # un-muting beta chain
            tool.mute_aa(0, 100000, "D")  # muting the alpha chain
            # Beta
            if verbose:
                subprocess.run([rosetta_dir, "-s", file_name, "@BETA_flags"])  # Runs InterfaceAnalyzer
            else:
                subprocess.run([rosetta_dir, "-s", file_name, "@BETA_flags"],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # Runs InterfaceAnalyzer
            tool.unmute_aa(0, 100000, "D")  # un-muting alpha chain
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


def write_results(results):
    """
    Writes the results of interface breakdown

    Parameters
    ----------
    results : dict
        information about the analysis of the interface scores of each chain and tcr
    """
    global verbose
    with open("AB_Usage.csv", "w") as f:
        f.write("PDB,ALPHA,BETA\n")
        if verbose:
            print("PDB\tALPHA\tBETA")
        for pdb_id in results:
            f.write(pdb_id + "," + str(results[pdb_id]["ALPHA"]) + "," + str(results[pdb_id]["BETA"]) + "\n")
            if verbose:  # Writes to stdout if verbose
                print(pdb_id + "\t" + str(results[pdb_id]["ALPHA"]) + "\t" + str(results[pdb_id]["BETA"]))


def write_flag(chains):
    """
    Creates flag files for InterfaceAnalyzer

    Parameters
    ----------
    chains : str
        chain of the TCR being analyzed
    """
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
# TODO: Option to save plot if not using GUI
def single_graph(score_file, x, y, ymax, ymin, xmax, xmin):
    """
    Produces a single seaborn relplot of the scoring of a native TCR redocking

    Parameters
    ----------
    score_file : str
        name and location of score file
    x : str
        value being plotted for x-axis
    y : str
        value being plotted for y-axi
    ymax : int
        max value for y axis
    ymin : int
        min value for y axi
    xmax : int
        max value for x axis
    xmin : int
        min value for x axis
    """
    content = subprocess.run("tr -s ' ' < " + score_file + " | tr ' ' ','", shell=True, stdout=subprocess.PIPE)
    if "/" in score_file:  # Detect if not in current dir
        dir_file = "/".join(score_file.split("/")[:-1]) + "/"
    else:
        dir_file = ""
    new_name = score_file.split("/")[-1].split(".")[0] + ".csv"
    with open(dir_file + new_name, "w") as f:
        first = True
        line_count = 0
        for line in content.stdout.decode('utf-8').split('\n')[1:]:
            if first:
                f.write(",".join(line.split(",")[1:]) + "\n")  # Remove SCORE: from first line
                first = False
            else:
                f.write(str(line_count) + "," + ",".join(line.split(",")[1:]) + '\n')
                line_count += 1
    score_df = pd.read_csv(dir_file + new_name, index_col=0)
    if ymax:
        score_df = score_df[score_df[y] <= ymax]
    if ymin:
        score_df = score_df[score_df[y] >= ymin]
    if xmax:
        score_df = score_df[score_df[x] <= xmax]
    if xmin:
        score_df = score_df[score_df[x] >= xmin]
    score_df = score_df.sort_values(y)
    colors = ["red"]  # Color lowest value
    colors += ["blue"] * (score_df["description"].count() - 1)
    sns.relplot(data=score_df, x=x, y=y, hue="description", palette=colors, legend=False)
    plt.show()
    os.remove(dir_file + new_name)


def multi_graph(dir_score_files, x, y, ymax, ymin, xmax, xmin):
    """
    Produces a multi-plot of seaborn relplots of the scoring of a native TCR redocking

    Parameters
    ----------
    dir_score_files : str
        directory that contains the .sc file
    x : str
        value being plotted for x-axis
    y : str
        value being plotted for y-axi
    ymax : int
        max value for y axis
    ymin : int
        min value for y axi
    xmax : int
        max value for x axis
    xmin : int
        min value for x axis
    """
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
        first = True
        line_count = 0
        for line in content.stdout.decode('utf-8').split('\n'):
            if first:
                f.write(",".join(line.split(",")[1:]) + "\n")  # Remove SCORE: from first line
                first = False
            else:
                f.write(line + '\n')
    score_df = pd.read_csv(dir_score_files + "/" + new_name, index_col=0)
    if ymax:
        score_df = score_df[score_df[y] <= ymax]
    if ymin:
        score_df = score_df[score_df[y] >= ymin]
    if xmax:
        score_df = score_df[score_df[x] <= xmax]
    if xmin:
        score_df = score_df[score_df[x] >= xmin]
    graph = sns.FacetGrid(score_df, col="file", hue="CAPRI_rank", col_wrap=4, palette="RdYlGn")
    # graph = sns.FacetGrid(score_df, col="file", hue="description", palette=colors, col_wrap=4)
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


def read_sheet(sheet_in):
    """
    Reading energy breakdown csv file
    File: 0: Score, 1: pose_id, 2: resi1, 3: pdbid1, 4: restype1, 5: resi2, 6: pdbid2, 7: restype2, 8: fa_atr
    9: fa_rep, 10: fa_sol, 11: fa_sol_rep, 12: fa_intra_sol_xover4, 13: ik_ball_wtd, 14: fa_elec, 15: pro_close
    16: hbond_sr_bb, 17: hbond_lr_bb, 18: hbond_bb_sc, 19: hbond_sc, 20: dslf_fa13, 21: omega, 22: fa_dun
    23: p_aa_pp, 24: yhh_planarity, 25: ref, 26: rama_prepro, 27: total, 28: description

    Parameters
    ----------
    sheet_in : str
        name / location of csv file

    Returns
    -------
    chains : dict
        {'A':[['1A','GLY],['2A','SER'],...]}
    interactions : dict
        {'1A':['1A', '2A',0.0],...],'2A':[..]}
    """
    interactions = {}  # Dir of interactions in PDB, Key:
    chains = {}  # Dir of chains in file, Key: Chain_ID, Value: 0 - red_id, 1 - AA
    df = pd.read_csv(sheet_in)
    for key, value in df.iterrows():
        values = value.array
        if values[7] != "onebody":
            if values[3][-1] != values[6][-1]:  # Ensure not capturing internal interactions NEW
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


def get_peptide_inter(chains, interactions, chain_in):
    """
    Gathers each interaction and produces a dictionary with {'1A':'GLY',...,'105D':'GLY'}

    Parameters
    ----------
    chains : dict
        Dictionary produced by above method - {'A':[['1A','GLY],['2A','SER'],...]}
    interactions :dict
        Dictionary produced by above method - {'1A':['1A', '2A',0.0],...],'2A':[..]}
    chain_in : str
        Chain being compared to for interactions to the TCR

    Returns
    -------
    aa_inter : dict
        Dictionary containing the interactions to be documented in the heatmap or csv table
    """
    peptide = []
    aa_inter = {}
    for AA_Peptide in chains[chain_in]:
        peptide.append(AA_Peptide)
    for each in peptide:
        if each[0] in interactions.keys():
            if each[0] in aa_inter.keys():
                aa_inter[each[0]].append(interactions[each[0]])
            else:
                aa_inter[each[0]] = [interactions[each[0]]]
    return aa_inter


def make_peptide_table(tcr_file, aa_list, chain_list, csv_name, chain_in):
    """
    Creates the CSV table of results of the interactions

    Parameters
    ----------
    tcr_file : str
        PDB of tcr being analyzed
    aa_list : dict
        List of interactions provided by method above
    chain_list : list
        List of chains being compared
    csv_name : str
        Name of the outputted table being produced
    chain_in : str
        Chain that's being determined what interactions it has to the TCR chains
    """
    aa_info = {}  # Key: ChainID Value: AA as 3 letter
    tool = PdbTools3(tcr_file)
    cdr_pos = tool.pull_cdr()
    cdr_info = {"D": {"CDR1A": range(cdr_pos[0][0][1], cdr_pos[0][0][2]),
                      "CDR2A": range(cdr_pos[0][1][1], cdr_pos[0][1][2]),
                      "CDR2.5A": range(cdr_pos[0][2][1], cdr_pos[0][2][2]),
                      "CDR3A": range(cdr_pos[0][3][1], cdr_pos[0][3][2])},
                "E": {"CDR1B": range(cdr_pos[1][0][1], cdr_pos[1][0][2]),
                      "CDR2B": range(cdr_pos[1][1][1], cdr_pos[1][1][2]),
                      "CDR2.5B": range(cdr_pos[1][2][1], cdr_pos[1][2][2]),
                      "CDR3B": range(cdr_pos[1][3][1], cdr_pos[1][3][2])}}
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


def heatmap_info(tcr_file, aa_list, chain_list, chain_in):
    """
    Generates data for heatmap based on chain_in and CDR regions
    Return csv with dataset to create heatmap & cdr_info to avoid rerunning search

    Parameters
    ----------
    tcr_file : str
        PDB of tcr being analyzed
    aa_list : dict
        List of interactions provided by method above
    chain_list : list
        List of chains being compared
    chain_in : str
        Chain that's being determined what interactions it has to the TCR chains

    Returns
    -------
    file_name : str
    cdr_info : dict
    """
    aa_info = {}  # Key: ChainID Value: aa as 3 letters
    tool = PdbTools3(tcr_file)
    cdr_pos = tool.pull_cdr()
    cdr_info = {"D": {"CDR1α": range(cdr_pos[0][0][1], cdr_pos[0][0][2]),
                      "CDR2α": range(cdr_pos[0][1][1], cdr_pos[0][1][2]),
                      "CDR2.5α": range(cdr_pos[0][2][1], cdr_pos[0][2][2]),
                      "CDR3α": range(cdr_pos[0][3][1], cdr_pos[0][3][2])},
                "E": {"CDR1β": range(cdr_pos[1][0][1], cdr_pos[1][0][2]),
                      "CDR2β": range(cdr_pos[1][1][1], cdr_pos[1][1][2]),
                      "CDR2.5β": range(cdr_pos[1][2][1], cdr_pos[1][2][2]),
                      "CDR3β": range(cdr_pos[1][3][1], cdr_pos[1][3][2])}}
    header_list = []  # Saves the header constructed: numAA
    for chain in chain_list:  # Key: '234E': 'VAL' for every AA
        for each in chain_list[chain]:
            aa_info[each[0]] = each[1]
    file_name = "temp.csv"
    with open(file_name, "w") as t2:
        for chain in cdr_info:  # First line: CDR id's
            for cdr in cdr_info[chain]:
                for num in cdr_info[chain][cdr]:
                    if (str(num) + chain) in aa_info:  # Catch when gap aa in range of CDR
                        header_list.append(str(num) + chain)
                        t2.write("," + aa_info[str(num) + chain] + " " + str(num))  # Prints header ex. SER 32, VAL 50..
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
    return file_name, cdr_info


def heatmap(info, cdr_info, chain_in, font_size, distance_in=0.0, distance=False):
    """
    Produce the heatmap with the informatino collected above using seaborn

    Parameters
    ----------
    info : str
        Location of CSV for heatmap creation
    cdr_info : dict
        Dictionary produced by above method
    chain_in : str
        Chain that's being determined what interactions it has to the TCR chains
    font_size : int
        Size the user wants the font to be
    distance_in : float
        vmax
    distance : boolean
        If utilizing distance cut on plot
    """
    # Read in csv
    df = pd.read_csv(info)
    # Remove columns with only zero values if MHC
    if chain_in == "A":
        # df = df.loc[:, (df != 0).any(axis=0)]
        df = df[(df.sum(axis=1) != 0)]
    # Read in x and y axis labels
    y_axis_labels = list(df.iloc[:,0])
    x_axis_labels = list(df.iloc[0:, :])[1:]
    # Drop first column
    df = df.iloc[:, 1:]
    # Annotation labels - only when below 0
    if not distance:
        labels = df.iloc[:, :].applymap(lambda v: str(v) if v < 0 else '')
        ax = sns.heatmap(df, vmax=0, linewidths=.2, linecolor="grey", xticklabels=x_axis_labels,
                         yticklabels=y_axis_labels,
                         cmap=sns.cubehelix_palette(start=2, rot=0, reverse=True, dark=0, light=1, as_cmap=True),
                         annot=labels, annot_kws={"fontsize": font_size, 'rotation': 90}, fmt='')
        ax.set_yticklabels(y_axis_labels, size=font_size)
        ax.set_xticklabels(x_axis_labels, size=font_size)
    if distance:
        labels = df.iloc[:, :].applymap(lambda v: str(f'{v:.3f}') if v < 10000 else '')
        ax = sns.heatmap(df, vmax=distance_in, vmin=0, linewidths=.2, linecolor="grey", xticklabels=x_axis_labels,
                         yticklabels=y_axis_labels,
                         cmap=(sns.cubehelix_palette(start=2, rot=0, reverse=True, dark=0, light=1)),
                         annot=labels, annot_kws={"fontsize": font_size, 'rotation': 90}, fmt='')
        ax.set_yticklabels(y_axis_labels, size=font_size)
        ax.set_xticklabels(x_axis_labels, size=font_size)
    plt.xlabel("TCR", size=font_size)
    y_label = "Peptide"
    if chain_in == "A":
        y_label = "MHC"
    plt.ylabel(y_label, size=font_size)
    # Adding in additional tick marks to account for labeling CDR regions
    # if not remove_0:
    ax2 = ax.twiny()
    ax2.set_xlim([0, ax2.get_xlim()[1]])
    ax2.set_xticks(ax.get_xticks())
    x_tick_labels = []
    # Capture first int
    chain_pos = ["D", "E"]
    chain_pos_num = 0
    cdr_pos = 0
    cdrs = ["CDR1α", "CDR2α", "CDR2.5α", "CDR3α", "CDR1β", "CDR2β", "CDR2.5β", "CDR3β"]
    cdr_color_dic = {"CDR1α": "red", "CDR2α": "orange", "CDR2.5α": "darkgreen", "CDR3α": "blue",
                     "CDR1β": "red", "CDR2β": "orange", "CDR2.5β": "darkgreen", "CDR3β": "blue"}
    colors = []
    results = []
    for each in x_axis_labels:
        num = int(each.split(" ")[-1].split(".")[0])
        if num not in cdr_info[chain_pos[chain_pos_num]][cdrs[cdr_pos]]:  # Choosing correct cdr color
            if cdr_pos == 3:  # Once you hit CDR3a go to CDR1b
                chain_pos_num = 1
            cdr_pos += 1  # Progress
        colors.append(cdr_color_dic[cdrs[cdr_pos]])
        results.append(cdrs[cdr_pos])
        x_tick_labels.append(cdrs[cdr_pos])
    # Set CDR xtick labels
    ax2.set_xticklabels(x_tick_labels, fontsize=font_size, rotation=90)
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


def interface_heatmap(tcr_file, interface_breakdown, mhc, font_size):
    """
    Controls the helper method to produce heatmap of TCR residue energy breakdown results

    Parameters
    ----------
    tcr_file : str
        Location of tcr file
    interface_breakdown : str
        Output of interface breakdown
    mhc : boolean
        If the user wants to compare the TCR interactions to the MHC chain
    font_size : int
        Size the user wants the font to be
    """
    # Changes to MHC chain if requested
    chain_in = "C"
    if mhc:
        chain_in = "A"
    chain_list, inter_list = read_sheet(interface_breakdown)
    aa_inter = get_peptide_inter(chain_list, inter_list, chain_in)
    info, cdr_info = heatmap_info(tcr_file, aa_inter, chain_list, chain_in)
    heatmap(info, cdr_info, chain_in, font_size)
    os.remove(info)


def get_contacts(tcr_file, distance, alpha_carbon, chain_in):
    """
    Calculate distances between CDR loops and antigen - save atoms within cutoff distance

    Parameters
    ----------
    tcr_file : str
        PDB of TCR
    distance : float
        Cutoff of collected values
    alpha_carbon : boolean
        Distance from alpha_carbon - else distance from closest atom in aa
    chain_in : chain_in
        MHC chain: A or Peptide chain: C

    Returns
    -------
    contacts : dict
        {AA comp_num: {partner AA comp_num: distance}}
    cdr_info : dict
        Dictionary containing information about CDR regions of each chains
    cdr_aa : dict
        Dictionary containing information about CDR regions for each chain
    all_aa : dict
        Store list of residues and resi num for each amino acid in each chain
    """
    tool = PdbTools3(tcr_file)
    cdr_pos = tool.pull_cdr()
    cdr_info = {"D": {"CDR1α": range(cdr_pos[0][0][1], cdr_pos[0][0][2]),
                      "CDR2α": range(cdr_pos[0][1][1], cdr_pos[0][1][2]),
                      "CDR2.5α": range(cdr_pos[0][2][1], cdr_pos[0][2][2]),
                      "CDR3α": range(cdr_pos[0][3][1], cdr_pos[0][3][2])},
                "E": {"CDR1β": range(cdr_pos[1][0][1], cdr_pos[1][0][2]),
                      "CDR2β": range(cdr_pos[1][1][1], cdr_pos[1][1][2]),
                      "CDR2.5β": range(cdr_pos[1][2][1], cdr_pos[1][2][2]),
                      "CDR3β": range(cdr_pos[1][3][1], cdr_pos[1][3][2])}}
    cdr_aa = {"D": [], "E": []}
    # Convert cdr_info into list of positions
    for chain in cdr_info:
        for cdr in cdr_info[chain]:
            cdr_aa[chain] += list(cdr_info[chain][cdr])
    antigen_atoms = tool.get_atoms_on_chain(chain_in)  # Collect antigen atoms based on chain_in
    if chain_in == "A":  # If pmhc, only collect first 180 amino acids
        aa_count = 0
        past_aa = -1
        temp_list = []
        for atom in antigen_atoms:
            if aa_count <= 180:
                if atom["comp_num"] != past_aa:
                    past_aa = atom["comp_num"]
                    aa_count += 1
                temp_list.append(atom)
        antigen_atoms = temp_list
    alpha_atoms = tool.get_atoms_on_chain("D")  # Collect alpha chain atoms
    beta_atoms = tool.get_atoms_on_chain("E")  # Collect beta chain atoms
    cdr_atoms = {"D": [], "E": []}  # Store list of atoms on CDR residues
    all_aa = {"D": {}, "E": {}, chain_in: {}}  # Store list of residues and resi num for each amino acid in each chain
    # Loop through all atoms
    for chain in [alpha_atoms, beta_atoms]:
        for atom in chain:
            if atom["comp_num"] not in list(all_aa[atom["chain_id"]]):  # Store all_aa info {chain: {202: THR}}
                all_aa[atom["chain_id"]][atom["comp_num"]] = atom["atom_comp_id"]
            if atom["comp_num"] in cdr_aa[atom["chain_id"]]:  # Find atoms in cdr loop
                cdr_atoms[atom["chain_id"]].append(atom)
    # print([atom for atom in cdr_atoms["D"] if atom["atom_id"]=="CA"])
    if alpha_carbon:
        for chain in cdr_atoms:
            cdr_atoms[chain] = [atom for atom in cdr_atoms[chain] if atom["atom_id"]=="CA"]
        antigen_atoms = [atom for atom in antigen_atoms if atom["atom_id"]=="CA"]
    contacts = {"D": {}, "E": {}}  # chain_id: {AA comp_num: {partner AA comp_num: distance}}
    for chain in cdr_atoms.keys():
        flag = False  # Once we've looped once
        for atom_1 in cdr_atoms[chain]:
            for atom_2 in antigen_atoms:
                # Collect antigen resi as above
                if not flag:
                    if atom_2["comp_num"] not in all_aa[atom_2["chain_id"]]:
                        all_aa[atom_2["chain_id"]][atom_2["comp_num"]] = atom_2["atom_comp_id"]
                # Calculate distance between each partner
                euc_dist = tool.euclidean_of_atoms(atom_1["atom_num"], atom_2["atom_num"])
                if euc_dist <= distance:  # if within threshold
                    if atom_1['comp_num'] in contacts[chain].keys():  # if cdr aa already documented
                        if atom_2['comp_num'] in contacts[chain][atom_1['comp_num']]:  # If cdr aa partner found
                            if euc_dist < contacts[chain][atom_1['comp_num']][atom_2['comp_num']]:
                                contacts[chain][atom_1['comp_num']][atom_2['comp_num']] = euc_dist
                        else:
                            contacts[chain][atom_1['comp_num']][atom_2['comp_num']] = euc_dist
                    else:  # if not previous aa to aa contact documented
                        contacts[chain][atom_1['comp_num']] = {atom_2['comp_num']: euc_dist}
            flag = True
    return contacts, cdr_info, cdr_aa, all_aa  # Chain_id: {aa comp_num: {partner AA comp_num: distance}}


def distance_heatmap_info(contacts, cdr_info, all_aa, chain_in):
    """
    Produce the heatmap csv information table

    Parameters
    ----------
    contacts : dict
        {AA comp_num: {partner AA comp_num: distance}}
    cdr_info : dict
        Dictionary containing information about CDR regions of each chains
    all_aa : dict
        Store list of residues and resi num for each amino acid in each chain
    chain_in : str

    Returns
    -------
    file_name : str
        Name of csv produced
    """
    file_name = "temp.csv"
    with open(file_name, "w") as f1:
        header_list = []  # Keeps tracks of position of aa in header list ex. 101E
        for chain in cdr_info:  # Alpha & beta chain
            for cdr in cdr_info[chain]:  # Loop over each cdr
                for num in cdr_info[chain][cdr]:  # Loop over each position in cdr
                    # Write Three Letter + Resi Num
                    if num in all_aa[chain].keys():
                        header_list.append(str(num) + chain)
                        f1.write("," + all_aa[chain][num] + " " + str(num))
        f1.write("\n")
        for num in all_aa[chain_in]:
            f1.write(all_aa[chain_in][num] + " " + str(num))
            for aa in header_list:
                if int(aa[:-1]) in contacts[aa[-1]].keys():
                    if num in contacts[aa[-1]][int(aa[:-1])].keys():  # Contacts[chain][num]
                        f1.write("," + str(contacts[aa[-1]][int(aa[:-1])][num]))
                    else:
                        f1.write("," + str(10000))
                else:
                    f1.write("," + str(10000))
            f1.write("\n")
    return file_name


def interface_heatmap_dist(tcr_file, distance, alpha_carbon, mhc, font_size):
    """
    Controls the methods used for interface heatmap with distance cutoff

    Parameters
    ----------
    tcr_file : str
        Location of TCR PDB file
    distance : float
        Distance cutoff
    alpha_carbon : boolean
        If using alpha_carbon distance
    mhc : boolean
        If user wants to compare TCR chains with MHC instead of peptide
    font_size : int
    """
    chain_in = "C"
    if mhc:
        chain_in = "A"
    contacts, cdr_info, cdr_aa, all_aa = get_contacts(tcr_file, distance, alpha_carbon, chain_in)  # Collect contacts
    temp_file = distance_heatmap_info(contacts, cdr_info, all_aa, chain_in)  # Create csv for heatmap
    heatmap(temp_file, cdr_info, chain_in, font_size, distance, True)  # Create heatmap
    os.remove(temp_file)


def peptide_table(tcr_file, energy_breakdown, mhc):
    """
    Controller method to generate energy breakdown table

    Parameters
    ----------
    tcr_file : str
        Location of TCR PDB file
    energy_breakdown : str
        Location fo energry breakdown output
    mhc : boolean
        If user wants to compare TCR chains with MHC instead of peptide
    """
    # Changes to MHC chain if requested
    chain_in = "C"
    if mhc:
        chain_in = "A"
    chain_list, inter_list = read_sheet(energy_breakdown)
    aa_inter = get_peptide_inter(chain_list, inter_list, chain_in)
    make_peptide_table(tcr_file, aa_inter, chain_list, "output.csv", chain_in)


def tsv_to_csv(tsv_in):
    """
    Convert tsv to csv

    Parameters
    ----------
    tsv_in : str
        Location of tsv file

    Returns
    -------
    new_name : str
        Updated location of file, now csv
    """
    name_in = tsv_in
    # Rename to same name but replace file extension to .csv
    new_name = "".join(tsv_in.split(".")[:-1]) + ".csv"
    submit_name = name_in.replace(" ", "\\ ")
    # Convert .out (tsv) to .csv
    cmd = "tr -s ' ' < " + submit_name + " | tr ' ' ','"
    new_file = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE)
    os.rename(name_in, new_name)
    with open(new_name, "w") as f1:
        f1.write(new_file.stdout.decode('utf-8')[1:])
    return new_name


def run_breakdown(pdb_in):
    """
    Run Rosetta Residue Energy Breakdown program if PDB submitted for heatmap or energy breakdown table

    Parameters
    ----------
    pdb_in : str
        Location of TCR PDB file

    Returns
    -------
    convert_out : str
        Location of the outputted csv results
    """
    # Update program location information
    global rosetta_dir
    # In and out files
    in_file = "-in:file:s " + pdb_in
    file_name = "energy_breakdown.out"
    out_file = "-out:file:silent " + file_name
    # Run process
    subprocess.run([rosetta_dir, in_file, out_file], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # Convert tsv to csv
    file_name = os.getcwd() + "/" + file_name
    convert_out = tsv_to_csv(file_name)
    return convert_out


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
    with open("../config.ini", "r") as f1:
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
            rosetta_dir += "main/source/bin/" + program
            break
        elif not program.startswith(program_in + ".default"):
            rosetta_dir += "main/source/bin/" + program


####################
#     Controls     #
####################
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--ab", help="(AB usage) Submit PDB for ab_usage interface scores")
    parser.add_argument("-v", "--verbose", help="(AB usage) Verbose", action="store_true", default=False)
    parser.add_argument("-s", "--sc", help="(Native) Score file produced from docking or refinement (or dir of .sc)")
    parser.add_argument("-x", "--xaxis", help="(Native) X axis for native structure comparison", type=str)
    parser.add_argument("-y", "--yaxis", help="(Native) Y axis for native structure comparison", type=str)
    parser.add_argument("--ymax", help="(Native) Y axis value maximum", type=int)
    parser.add_argument("--ymin", help="(Native) Y axis value minimum", type=int)
    parser.add_argument("--xmax", help="(Native) X axis value maximum", type=int)
    parser.add_argument("--xmin", help="(Native) X axis value minimum", type=int)
    parser.add_argument("-p", "--heatmap_dist", help="(DHM) Distance Breakdown TCR PDB File", type=str)
    parser.add_argument("-d", "--distance", help="(DHM) Distance for heatmap | Default 4.5", type=float, default=4.5)
    parser.add_argument("-c", "--alpha_carbon", help="(DHM) Alpha Carbon Only", action="store_true", default=False)
    parser.add_argument("-e", "--heatmap_energy", help="(EB) TCR PDB File", type=str)
    parser.add_argument("-t", "--table", help="(EB) TCR PDB File", type=str)
    parser.add_argument("-m", "--mhc", help="(DHM/EB) Changes energy breakdown to MHC versus peptide",
                        action="store_true", default=False)
    parser.add_argument("-f", "--fontsize", help="(DHM/EB) Adjust font size | Default 12", type=int, default=12)
    return parser.parse_args()


def main():
    args = parse_args()
    # AB usage output
    if args.ab:
        global verbose
        if args.verbose:
            global verbose
            verbose = True
        rosetta_binary("InterfaceAnalyzer.")  # Set binary
        # Set TCR chains based on user input
        tcr_chains = {"MHC": "A", "peptide": "C",
                      "ALPHA": "D", "BETA": "E"}
        run_interface(args.ab, tcr_chains)
    # Native structure comparison
    if args.sc:
        if not os.path.isdir(args.sc):
            single_graph(args.sc, args.xaxis, args.yaxis, args.ymax, args.ymin, args.xmax, args.xmin)
        else:
            multi_graph(args.sc, args.xaxis, args.yaxis, args.ymax, args.ymin, args.xmax, args.xmin)
    # Heatmap or table routing for energy breakdown
    if args.heatmap_energy:  # Generate heatmap
        if args.heatmap_energy.endswith(".pdb"):  # Run energy breakdown if PDB submitted
            rosetta_binary("residue_energy_breakdown")  # Set binary
            breakdown_file = run_breakdown(args.heatmap_energy)
            interface_heatmap(args.heatmap_energy, breakdown_file, args.mhc, args.fontsize)
            os.remove(breakdown_file)
    if args.table:  # Generate breakdown table
        if args.table.endswith(".pdb"):  # Run energy breakdown if PDB submitted
            rosetta_binary("residue_energy_breakdown")  # Set binary
            breakdown_file = run_breakdown(args.table)
            peptide_table(args.table, breakdown_file, args.mhc)
            os.remove(breakdown_file)
    # Heatmap for distance
    if args.heatmap_dist:
        if args.heatmap_dist.endswith(".pdb"):
            interface_heatmap_dist(args.heatmap_dist, args.distance, args.alpha_carbon, args.mhc, args.fontsize)
        else:
            print("Provide PDB File")


if __name__ == '__main__':
    main()
