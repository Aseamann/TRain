# This file is a part of the TRain program
# Author: Austin Seamann & Dario Ghersi
# Version: 0.30
# Last Updated: October 22, 2021
import argparse
import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import subprocess


#################
#    Methods    #
#################
def density(df_0):
    df = pd.read_csv(df_0, index_col=0, sep="\t")
    ax_1 = sns.displot(df, x="I_sc", hue="TCR", kind="kde", fill=True)
    ax_1.set(title="TCR I_sc Density")
    df_2 = []
    for tcr in df["TCR"]:
        df_2.append(df[(df.tcr == tcr) and (df.pmhc == tcr)])
    plt.plot(df_2)
    plt.show()
    # for tcr in df["TCR"].unique():
    #     # df_1 = df.drop(df.index[df['TCR'] != tcr], inplace=False)
    #     # df["I_sc"] = -1 * df["I_sc"]
    #     # ax = sns.displot(df_1, x="I_sc", hue="pMHC", element="step", binwidth=2)
    #     ax.set(title=tcr)
    #     plt.show()


def strip_plot_isc(df_0):
    score_df = pd.read_csv(df_0, index_col=0, sep="\t").sort_values('TCR')
    ebv_df = score_df.query("pMHC_Ant == 'EBV'")
    m1_df = score_df.query("pMHC_Ant == 'M1'")
    tax_df = score_df.query("pMHC_Ant == 'Tax'")
    sns.stripplot(x="TCR", y="I_sc", data=score_df, hue="pMHC_Ant").set(title="TCR I_sc Strip Plot")
    sns.lineplot(x="TCR", y="I_sc", data=m1_df)
    sns.lineplot(x="TCR", y="I_sc", data=ebv_df)
    sns.lineplot(x="TCR", y="I_sc", data=tax_df)
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
    plt.show()


def strip_plot_ab(df_0):
    score_df = pd.read_csv(df_0, index_col=0, sep="\t").sort_values("TCR")
    temp_ab = []
    for index, row in score_df.iterrows():
        temp_ab.append(float(row['alpha']) + float(row['beta']))
    score_df['total'] = temp_ab
    print(score_df)
    sns.stripplot(x="TCR", y="total", data=score_df, hue="pMHC_Ant").set(title="TCR Total Score Strip Plot")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
    plt.show()


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


def i_v_rmsd(score_file):
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
    sns.relplot(data=score_df, x="rms", y="I_sc")
    plt.show()


def t_v_rmsd(score_file):
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
    sns.relplot(data=score_df, x="rms", y="total_score")
    plt.show()


def mds_v_rmsd(score_file):
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
    sns.relplot(data=score_df, x="rms", y="motif_dock")
    plt.show()


def multi_x_rmsd(dir_score_files, x):
    compare = {"irmsd": "I_sc", "mdsrmsd": "motif_dock", "trmsd": "total_score"}
    all_file = os.getcwd() + "/" + dir_score_files + "/all.sc"
    print(all_file)
    header = False
    current_file = ""
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
    graph = sns.FacetGrid(score_df, col="file", hue="CAPRI_rank", col_wrap=4, palette="RdYlGn", xlim=(0, 200))
    graph.map(sns.scatterplot, "cen_rms", compare[x])
    graph.add_legend()
    plt.show()


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
    with open(csv_name, "w") as t1:
        t1.write("Chain: " + chain_in + ",Interaction energy,Partner,CDR Region\n")
        for AA in aa_list:
            for each in aa_list[AA][0]:
                if each[0] == AA:
                    partner = each[1]
                else:
                    partner = each[0]
                if partner[-1] == "D" or partner[-1] == "E":  # Only allows for TCR chains
                    output = aa_info[AA] + " " + str(int(AA[:-1]) - int(first_aa[AA[-1]]) + 1) + "," + str(each[2]) \
                             + "," + aa_info[partner] + " " + partner[:-1] + "," + partner[-1] + "\n"
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


def heatmap(info, remove_0):
    # Read in csv
    df = pd.read_csv(info)
    # Remove columns with only zero values
    if remove_0:
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
    plt.ylabel("Peptide")
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


def interface_heatmap(interface_breakdown):
    chain_list, inter_list = read_sheet(interface_breakdown)
    aa_inter = get_peptide_inter(chain_list, inter_list, "A")
    info = heatmap_info(aa_inter, chain_list, "A")
    heatmap(info, True)
    os.remove(info)


####################
#     Controls     #
####################
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", help="TSV containing: pdb_name,TCR,pMHC,Total_Score,Alpha,Beta", type=str)
    parser.add_argument("-s", "--strip", help="TSV containing: pdb_name,TCR,pMHC,Total_Score,Alpha,Beta")
    parser.add_argument("-z", "--strip2", help="TSV containing: pdb_name,TCR,pMHC,Total_Score,Alpha,Beta"\
                        " | produces strip plot based on alpha/beta score")
    parser.add_argument("-b", "--box_tcr", help="Box plot of TCR in comparison of native to all others.",
                        default=False, action="store_true")
    parser.add_argument("-i", "--irmsd", help="Produces RMSD vs I_sc for inputted .sc file", type=str)
    parser.add_argument("-m", "--mdsrmsd", help="Produces RMSD vs MDS for inputted .sc file", type=str)
    parser.add_argument("-t", "--trmsd", help="Produces RMSD vs Total_score for inputted .sc file", type=str)
    parser.add_argument("-a", "--allscore", help="True if wanting to compare a directory of .sc files",
                        action="store_true", default=False)
    parser.add_argument("--heatmap", help="Energy breakdown csv", type=str)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.strip:
        strip_plot_isc(args.strip)
    if args.strip2:
        strip_plot_ab(args.strip2)
    if args.box_tcr:
        box_plot_isc(args.data)
    if args.irmsd:
        if not args.allscore:
            i_v_rmsd(args.irmsd)
        else:
            multi_x_rmsd(args.irmsd, "irmsd")
    if args.mdsrmsd:
        if not args.allscore:
            mds_v_rmsd(args.mdsrmsd)
        else:
            multi_x_rmsd(args.mdsrmsd, "mdsrmsd")
    if args.trmsd:
        if not args.allscore:
            t_v_rmsd(args.trmsd)
        else:
            multi_x_rmsd(args.trmsd, "trmsd")
    if args.heatmap:
        interface_heatmap(args.heatmap)


if __name__ == '__main__':
    main()
