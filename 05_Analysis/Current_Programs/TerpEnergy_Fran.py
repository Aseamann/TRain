import pandas
import argparse
import sys


first_aa = {}


# Reading energy breakdown excel file
# File: 0: Score, 1: pose_id, 2: resi1, 3: pdbid1, 4: restype1, 5: resi2, 6: pdbid2, 7: restype2, 8: fa_atr
# 9: fa_rep, 10: fa_sol, 11: fa_sol_rep, 12: fa_intra_sol_xover4, 13: ik_ball_wtd, 14: fa_elec, 15: pro_close
# 16: hbond_sr_bb, 17: hbond_lr_bb, 18: hbond_bb_sc, 19: hbond_sc, 20: dslf_fa13, 21: omega, 22: fa_dun
# 23: p_aa_pp, 24: yhh_planarity, 25: ref, 26: rama_prepro, 27: total, 28: description
def read_sheet(sheet_in):
    interactions = {}  # Dir of interactions in PDB, Key:
    chains = {}  # Dir of chains in file, Key: Chain_ID, Value: 0 - red_id, 1 - AA
    df = pandas.read_csv(sheet_in)
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


def get_peptide_inter(chains, interactions):
    peptide = []
    aa_inter = {}
    for AA_Peptide in chains["C"]:
        peptide.append(AA_Peptide)
    for each in peptide:
        if each[0] in aa_inter.keys():
            aa_inter[each[0]].append(interactions[each[0]])
        else:
            aa_inter[each[0]] = [interactions[each[0]]]
    return aa_inter


# Creates the CSV table
def make_peptide_table(aa_list, chain_list, csv_name):
    aa_info = {}  # Key: ChainID Value: AA as 3 letter
    cdr_info = {"D": {"CDR1A": range(28, 33), "CDR2A": range(50, 56), "CDR3A": range(91, 100)},
                "E": {"CDR1B": range(29, 33), "CDR2B": range(51, 59), "CDR3B": range(94, 104)}}
    for chain in chain_list:
        for each in chain_list[chain]:
            aa_info[each[0]] = each[1]
    with open(csv_name, "w") as t1:
        t1.write("Peptide,Interaction energy,Partner,CDR Region\n")
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
                        # print(cdr_info[partner[-1]][cdr])
                        if int(partner[:-1]) in cdr_info[partner[-1]][cdr]:
                            output = output[:-3] + "," + cdr + "\n"
                    t1.write(output)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Inputted energy xlsx file", type=str)
    parser.add_argument("-o", "--output", help="Output file name", default=sys.stdout, type=str)
    return parser.parse_args()


def main():
    args = parse_args()
    chain_list, inter_list = read_sheet(args.input)
    aa_inter = get_peptide_inter(chain_list, inter_list)
    # print(chain_list)
    # print(inter_list)
    make_peptide_table(aa_inter, chain_list, args.output)


if __name__ == "__main__":
    main()

