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
def get_tcr_info(file, sheet, positions):
    tcr_dic = {}  # Creates output dictionary
    clone_id_count = {}  # Keeps track of repeat clone ids
    previous_id = ""
    # if not clone_id provided
    if len(positions) == 6:
        positions.insert(0, -1)
    if file.endswith(".xlsx"):  # Open if xlsx
        if sheet == "0":  # If no sheet name was provided, set default to primary sheet
            sheet = 0
        df = pandas.read_excel(file, sheet_name=sheet)
    elif file.endswith(".csv"):  # Open if csv
        df = pandas.read_csv(file)
    for key, value in df.iterrows():  # Reading over table with pandas
        values = value.array
        # Finds each clone_id.
        if isinstance(values[positions[1]], str):
            # If no clone id, update with previous name or no clone id provided
            if not isinstance(values[positions[0]], str) or positions[0] == -1:
                if previous_id != "":
                    values[positions[0]] = previous_id
                else:  # When not using clone id
                    value[positions[0]] = "id-000"
            # Remove second class (ex.  TRAV12-2*01-*03 and TRBV20-1*01-*07: TRAV12-2*01 and TRBV20-1*01)
            # .split(";")[0] splits for when record have two gene segments listed, typically identical segments in aa
            if values[positions[1]].count("-") >= 2:
                values[positions[1]] = "-".join(values[positions[1]].split("-")[:-1])
            if values[positions[4]].split(";")[0].count("-") >= 2:
                values[positions[4]] = "-".join(values[positions[4]].split(";")[0].split("-")[:-1])
            # Confirm each gene segment name is in-fact a name in gene_dic
            if confirm_segment([values[positions[1]], values[positions[3]], values[positions[4]].split(";")[0],
                                values[positions[6]]]):
                if values[positions[0]] in clone_id_count.keys():  # Catch duplicate clone id names and rename "-01"
                    clone_id_count[values[positions[0]]] += 1
                    values[positions[0]] = \
                        f"{values[positions[0]]}-{str(clone_id_count[values[positions[0]]]).zfill(3)}"
                else:
                    clone_id_count[values[positions[0]].split("-")[0]] = 0
                previous_id = values[positions[0]].split("-")[0]
                tcr_dic[values[positions[0]]] = [values[positions[1]], values[positions[2]], values[positions[3]],
                                                 values[positions[4]].split(";")[0], values[positions[5]],
                                                 values[positions[6]]]
            else:
                print(f"Skipped due to gene segment not found:\n{value}")
    return tcr_dic


# Method: 03 - confirm_segment()
# Input: List of gene segments w/ clone id position 0
# Goal: Loop through gene segments of a clone id and return if a segment is not found in gene_dic
# Output: Clone ID if gene segment not discovered
def confirm_segment(seg_list):
    for seg in seg_list:
        if seg not in gene_dic.keys():
            return False
    return True


# Method: 04 - make_tcr_seq()
# Input: tcr_parts - dictionary containing each aa segment of TCR as constructed from get_tcr_info
# Goal: Sends tcr parts into align_overlap to create resulting Alpha and Beta chain sequences.
# Output: a dictionary with clone ID key and a list of Alpha and Beta chain sequences.
def make_tcr_seq(tcr_parts):
    output_seqs = {}
    # For each clone-id
    for each in tcr_parts:
        # Alignment construction for TRAV, and CDR3A
        temp_seq_alpha = align_overlap(gene_dic[tcr_parts[each][0]], tcr_parts[each][1], "V")
        # Alignment construction for TRAV + CDR3A and TRAJ
        temp_seq_alpha = str(align_overlap(temp_seq_alpha, gene_dic[tcr_parts[each][2]], "J"))
        # # Alignment construction for TRBV, and CDR3B
        temp_seq_beta = align_overlap(gene_dic[tcr_parts[each][3]], tcr_parts[each][4], "V")
        # Alignment construction for TRBV + CDR3B and TRBJ
        temp_seq_beta = align_overlap(temp_seq_beta, gene_dic[tcr_parts[each][5]], "J")
        parts = tcr_parts[each]
        output_seqs[each] = [temp_seq_alpha, temp_seq_beta, parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]]
    return output_seqs


# Method: 05 - align_overlap()
# Input:
#   front - start of seq
#   end - end of seq
#   part - If overlapping V or J segment
# Goal: Looks for overlap between front and end AAs besides when single AA overlap
# Output: Resulting overlapped segments -- if == J: full TCR sequence
def align_overlap(front, end, part):
    matrix = matlist.blosum62  # Protein alignment scoring matrix
    score = -1
    best_align = []
    if part == "V":
        overlap_cut = 4
    else:
        overlap_cut = 0
    # Finding shortest alignment with best score to append regions
    # overlap cut modified for v or j based on testing
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
        best_align = [temp_align[0][0], temp_align[0][1], temp_align[0][4]]  # Start seq, end seq, end position
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


# Method: 06 - if_gap()
# Input:
#   front -
#   end -
#   region -
#   gap_count -
# Goal:
# Output: Find alignment of segments
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


# Method: 07 - append_constant
# Input: tcr_seq_dic - dictionary of tcr components clone_id: alpha, beta, [other parts]
# Output: Updated tcr_seq_dic with alpha and beta chains replaced with concatenated constant regions
def append_constant(tcr_seq_dic, organism):
    constant_dic = create_constant_dic(organism)
    for tcr_id in tcr_seq_dic:
        tcr_seq_dic[tcr_id][0] == "ALPHA"
        tcr_seq_dic[tcr_id][1] == "BETA"


# Method: 08 - create_gene_dic()
# Input: fasta file containing IMGT data of constant region gene segments of TCRs
# Goal: Populate the gene dictionary from IMGT data.
# Output:
#   constant_dic =
def create_constant_dic(organism):
    constant_dic = {}
    with open("constant_regions.fasta", "r") as f:
        temp_gene_segment = ""
        for line in f:
            if line.split("|")[2] == organism:
                if line[0] == ">":  # Find header
                    header = line.split("|")
                    temp_imgt_id = header[0][1:]
                    constant_dic[temp_imgt_id] = []  # prepares location in dictionary
                # Collecting sequence data from fasta file per header
                else:
                    if line.__contains__('*'):
                        constant_dic[temp_fam] += line[:-1].replace('*', "")
                    elif line.__contains__('X'):
                        constant_dic[temp_fam] += line[:-1].replace('X', "")
                    else:
                        constant_dic[temp_fam] += line[:-1]
        else:
            print("Organism not found (Appending constant region)")
    return constant_dic


# Method 07: make_info_table()
# Input:
#   tcr_dic - resulting clone ids and tcr information
#   file_in - optional file name for information table
# Output: Table with clone id, AV, CDR3a, AJ, BV, CDR3b, BJ, full alpha seq, full beta seq (csv)
def make_info_table(tcr_dic, file_in="translated_seq.csv"):
    with open(file_in, "w") as f:  # Open file, name is set by default but can be changed
        f.write("clone ID,TRAV,CDR3 Alpha,TRAJ,TRBV,CDR3 Beta,TRBJ,ALPHA,BETA\n")  # Header
        for each in tcr_dic:  # Loop through each clone ID
            info = tcr_dic[each][2:]
            chains = tcr_dic[each][:2]
            f.write(each + "," + ",".join(info) + "," + ",".join(chains) + "\n")


# Method 08: make_rep_builder_file()
# Input:
#   alpha_file - name of alpha chain file
#   beta_file - name of beta chain file
#   tcr_dic - dictionary of all tcrs
# Output:
#   alpha_file - outputted fasta file for alpha chains for rep builder
#   beta_file - outputted fasta file for beta chains for rep builder
def make_rep_tcrmod_file(alpha_file, beta_file, tcr_dic):
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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("single_cell_table", help="Xlsx or Csv of single cell t-cell data", type=str)
    parser.add_argument("-s", "--sheet", help="Sheet name in table", default="0", type=str)
    parser.add_argument("-r", "--rep_builder", help="(Output) Repertoire Builder input file format", default=False,
                        action="store_true")
    parser.add_argument("-l", "--lyra", help="(Output) Lyra input file format", default=False, action="store_true")
    parser.add_argument("-t", "--tcrmodel", help="(Output) TCRmodel input file format", default=False,
                        action="store_true")
    parser.add_argument("-i", "--information", help="(Output) Creates information table of constructed TCR seq."
                        , default=False , action="store_true")
    parser.add_argument("--genefamily", help="Replacement gene family segment file", default="family_seq.fasta"
                        , type=str)
    parser.add_argument("-c", "--columns",
                        help="Positions of gene segments: Clone ID,AV,CDR3a,AJ,BV,CDR3b,BJ | Can exclude Clone ID",
                        default="3,4,5,7,8,9,11", type=str)
    parser.add_argument("-f", "--full", help="Append each chain with constant region | provide organism", type=str)
    return parser.parse_args()


def main():
    args = parse_args()  # collect arguments
    # create dictionary for program to pull gene family information from
    create_gene_dic(args.genefamily)
    # create tcr dictionary containing each segment
    tcr_dic = get_tcr_info(args.single_cell_table, args.sheet, [int(i) for i in args.columns.split(",")])
    # create full tcr sequences
    tcr_seq_dic = make_tcr_seq(tcr_dic)
    if args.full:
        tcr_seq_dic = append_constant(tcr_seq_dic, args.full)
    # Create information table with full TCR sequences
    if args.information:
        make_info_table(tcr_seq_dic, ".".join(args.single_cell_table.split(".")[:-1]) + "_Results.csv")
    if args.rep_builder or args.tcrmodel:
        make_rep_tcrmod_file("alpha.fasta", "beta.fasta", tcr_seq_dic)


if __name__ == '__main__':
    main()
