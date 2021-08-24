# This file is a part of the TRain program
# Author: Austin Seamann & Dario Ghersi
# Version: 0.02
# Last Updated: Aug 23rd, 2021
import argparse
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


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
    return parser.parse_args()


def main():
    args = parse_args()
    if args.strip:
        strip_plot_isc(args.strip)
    if args.strip2:
        strip_plot_ab(args.strip2)
    if args.box_tcr:
        box_plot_isc(args.data)


if __name__ == '__main__':
    main()
