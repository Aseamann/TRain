# This file is a part of the TRain program
# Author: Austin Seamann & Dario Ghersi
# Version: 1.0
# Last Updated: July 11th, 2021

import argparse
import os
import pandas
from Bio import pairwise2
import warnings
from Bio import BiopythonDeprecationWarning
with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
    from Bio.SubsMat import MatrixInfo as matlist


#################
# Global
#################
gene_dic = {}


#################
# Methods
#################
# Method: 01 - create_gene_dic()
# Input: fasta file containing IMGT data of gene segments of TCRs
# Goal: Populate the gene dictionary from IMGT data.
# Output:
#   gene_dic: Gene Segment Name: AA Sequence
def create_gene_dic(fasta_file):
    with open(fasta_file, "r") as family:
        temp_fam = ""  # Holds current gene segment name
        for line in family:
            if line[0] == ">":  # Find header
                header = line.split("|")
                temp_fam = header[1]  # gene segment name, creates temp header name
                gene_dic[temp_fam] = ""  # prepares location in dictionary
            # Collecting sequence data from fasta file per header
            else:
                if line.__contains__('*'):
                    gene_dic[temp_fam] += line[:-1].replace('*', "")
                else:
                    gene_dic[temp_fam] += line[:-1]
    return gene_dic


# Method: 02 - get_tcr_info()
# Input:
#   file: .xlsx or .csv file containing AV, AJ, BV, BJ, CDR3A, CDR3B
#   sheet: optional parameter for .xlsx file to collect data from a secondary sheet (default="Sheet1")
#   default positions: 3: cloneID, 4: AV, 5: CDR3A, 7: AJ, 8: BV, 9: CDR3B, 11: BJ
# Goal: Create dictionary of components of TCR chains from table
# Output:
#   tcr_dic: clone ID Positions: 0 - TRAV, 1 - CDR3A, 2 - TRAJ, 3 - TRBV, 4 - CDR3B, 5 - TRBJ, clone ID - 6
def get_tcr_info(file, sheet):
    tcr_dic = {}  # Creates output dictionary
    non_cdr = "LIST OF POSITIONS NOT CONTAINING CDRs"
    if file.endswith(".xlsx"):  # Open if xlsx
        df = pandas.read_excel(file, sheet_name=sheet)
    elif file.endswith(".csv"):  # Open if csv
        df = pandas.read_csv(file)
    for key, value in df.iterrows():  # Reading over table with pandas
        values = value.array
        # Catches right now for cloneIDs that are missing components.
        if isinstance(values[3], str):
            # Remove second class (ex.  TRAV12-2*01-*03 and TRBV20-1*01-*07: TRAV12-2*01 and TRBV20-1*01)
            if values[4].count("-") >= 2:
                values[4] = "-".join(values[4].split("-")[:-1])
            if values[8].count("-") >= 2:
                values[8] = "-".join(values[8].split("-")[:-1])
            # Catches for when VRBJ has two gene families listed for value 8
            tcr_dic[values[3]] = [values[4], values[5], values[7], values[8].split(";")[0], values[9], values[11]]
        else:
            print("Skipped due to gene segment not found: " + str(values[3]))
    return tcr_dic


# Sends tcr parts into align_overlap to create resulting Alpha and Beta chain sequences.
# Outputs a dictionary with clone ID key and a list of Alpha and Beta chain sequences.
def make_tcr_seq(tcr_parts):
    output_seqs = {}
    count = 0
    # For each clone-id
    for each in tcr_parts:
        print(each)
        print("ALPHA:")
        # Alignment construction for TRAV, and CDR3A
        temp_seq_alpha = align_overlap(gene_dic[tcr_parts[each][0]], tcr_parts[each][1], "V")
        print("Va + CDR3a:")
        print(temp_seq_alpha)
        # Alignment construction for TRAV + CDR3A and TRAJ
        print("Ja:")
        print(gene_dic[tcr_parts[each][2]])
        temp_seq_alpha = str(align_overlap(temp_seq_alpha, gene_dic[tcr_parts[each][2]], "J"))
        print("Result: ")
        print(temp_seq_alpha)
        # # Alignment construction for TRBV, and CDR3B
        temp_seq_beta = align_overlap(gene_dic[tcr_parts[each][3]], tcr_parts[each][4], "V")
        print("BETA:")
        print("Vb + CDR3b:")
        print(temp_seq_beta)
        print("Jb:")
        print(gene_dic[tcr_parts[each][5]])
        # Alignment construction for TRBV + CDR3B and TRBJ
        temp_seq_beta = align_overlap(temp_seq_beta, gene_dic[tcr_parts[each][5]], "J")
        print("Result:")
        print(temp_seq_beta)
        print("\n")
        parts = tcr_parts[each]
        output_seqs[each] = [temp_seq_alpha, temp_seq_beta, parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]]
        if count > 10:
            break
        else:
            count += 1
    return output_seqs


# Looks for overlap between front and end AAs besides when single AA overlap
def align_overlap(front, end, part):
    matrix = matlist.blosum62  # Protein alignment scoring matrix
    score = -1
    best_align = []
    if part == "V":
        overlap_cut = 4
    else:
        overlap_cut = 0
    # Finding shortest alignment with best score to append regions
    #overlap cut modified for v or j based on testing
    if part == "V":
        for length in range(len(end) - overlap_cut, 1, -1):
            temp_align = pairwise2.align.localdx(front[length * -1:], end, matrix, one_alignment_only=True)
            if temp_align[0][3] == 0:
                best_align = []
                score = temp_align[0][2]
                best_align.append(temp_align[0][0])  # Saves front chain
                best_align.append(temp_align[0][1])  # Saves end chain
                best_align.append(length)
    elif part == "J":
        offset = 6  # Helps with alignments to avoid too spaced out alignment
        temp_align = pairwise2.align.localms(front[(len(end) - offset) * -1:], end, 2, -1, -1, -0.5
                                             , one_alignment_only=True, penalize_extend_when_opening=False)
        print("Temp_align:")
        print(temp_align)
        best_align = [temp_align[0][0], temp_align[0][1], temp_align[0][4]]  # Start seq, end seq, end position
    print("Best align:")
    print(best_align)
    print("\n")
    if part == "V":
        if score == -1:  # When score is below 0
            output = str(front + end)
            return output
        elif best_align[0][:best_align[2]] == best_align[1][:best_align[2]]:  # When overlap is not trimming either seq.
            output = str(front[:best_align[2] * -1] + end)
            return output
        else:  # When score is above 0
            par_seq = if_gap(best_align[0], end, part, 0)
            output = front[:best_align[2] * -1] + par_seq
            return output
    elif part == "J":  # If joining region overlap
        output = front + best_align[1][best_align[2]:]
        return output


def count_end_gap(seq_in):
    count = 0
    flag = True  # For more than one gap
    for each in seq_in:
        if each == "-":
            count += 1
        else:
            count = 0
    return count


def if_break(seq_in):
    previous1 = ""
    previous2 = ""
    for letter in seq_in:
        if letter != "-":
            if previous1 == "-":
                if previous2 != "-":
                    return True
        previous2 = previous1
        previous1 = letter
    return False


def if_gap(front, end, region, gap_count):  # Check for cdr len in previous method, send in if V or J region.
    temp_seq = ""
    if region == "V":  # When the variable region and cdr are being overlapped.
        for letter in front:
            if letter != "-":
                temp_seq += letter
            else:
                return temp_seq + end[len(temp_seq):]
    elif region == "J":  # When the joining region and cdr are being overlapped.
        return front + end[gap_count:]


def make_rep_builder_file(alpha_file, beta_file, tcr_dic):
    with open(alpha_file, "w") as file_a:
        with open(beta_file, "w") as file_b:
            for clone_id in tcr_dic:
                tcr_alpha_chain = tcr_dic[clone_id][0]
                tcr_beta_chain = tcr_dic[clone_id][1]
                count_1 = 1
                count_2 = 1
                if tcr_alpha_chain != '':
                    file_a.write('>' + clone_id + "\n")
                    for aa in tcr_alpha_chain:
                        if count_1 % 81 != 0:
                            file_a.write(aa)
                            count_1 += 1
                        else:
                            file_a.write('\n' + aa)
                            count_1 = 2
                    file_a.write('\n')
                if tcr_beta_chain != '':
                    file_b.write('>' + clone_id + '\n')
                    for aa in tcr_beta_chain:
                        if count_2 % 81 != 0:
                            file_b.write(aa)
                            count_2 += 1
                        else:
                            file_b.write('\n' + aa)
                            count_2 = 2
                    file_b.write('\n')


def make_info_table(tcr_dic, file_in="translated_seq.csv"):
    with open(file_in, "w") as f:  # Open file, name is set by default but can be changed
        f.write("clone ID,TRAV,CDR3 Alpha,TRAJ,TRBV,CDR3 Beta,TRBJ,ALPHA,BETA\n")  # Header
        for each in tcr_dic:  # Loop through each clone ID
            info = tcr_dic[each][2:]
            chains = tcr_dic[each][:2]
            f.write(each + "," + ",".join(info) + "," + ",".join(chains) + "\n")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("single_cell_table", help="Xlsx or Csv of single cell t-cell data", type=str)
    parser.add_argument("-s", "--sheet", help="Sheet name in table", default="Sheet1", type=str)
    parser.add_argument("-r", "--rep_builder", help="(Output) Repertoire Builder input file format", default=""
                        , type=str)
    parser.add_argument("-l", "--lyra", help="(Output) Lyra input file format", default="", type=str)
    parser.add_argument("-t", "--tcrmodel", help="(Output) TCRmodel input file format", default="", type=str)
    parser.add_argument("-i", "--information", help="(Output) Creates information table of constructed TCR seq."
                        , default=False , action="store_true")
    parser.add_argument("--genefamily", help="Replacement gene family segment file", default="family_seq.fasta"
                        , type=str)
    return parser.parse_args()


def main():
    args = parse_args()  # collect arguments
    create_gene_dic(args.genefamily)  # create dictionary for program to pull gene family information from
    tcr_dic = get_tcr_info(args.single_cell_table, args.sheet)  # create tcr dictionary containing each segment
    tcr_seq_dic = make_tcr_seq(tcr_dic)  # create full tcr sequences
    if args.information:  # Create information table with full TCR sequences
        make_info_table(tcr_seq_dic, ".".join(args.single_cell_table.split(".")[:-1]) + ".csv")
    if args.rep_builder:
        make_rep_builder_file("alpha_test.fasta", "beta_test.fasta", tcr_seq_dic)


if __name__ == '__main__':
    main()
