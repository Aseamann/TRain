#!/usr/bin/python3

######################################################################
# SeqConductor.py -- A component of TRain                            #
# Copyright: Austin Seamann & Dario Ghersi                           #
# Version: 0.1                                                       #
# Last Updated: April 27th, 2022                                     #
# Goal: Take in gene segments and CDR3 regions of TCR chains and     #
#        produce full amino acid sequences                           #
#                                                                    #
# Positional argument: Single Cell Table (XLSX or CSV)               #
# Named arguments: -s --sheet SHEET NAME (Alternative sheet name if  #
#                             not primary sheet in XLSX file)        #
#                  -f --fasta (Output option -- Fasta output option) #
#                  -t --information (Output option -- CSV output     #
#                                   for information table results)   #
#                  -g --genefamily (Replacement gene family segment  #
#                                   if not homo sapiens)             #
#                  -c --columns COMMA SEP LIST (Positions of gene    #
#                               segments: Clone ID, AV, CDR3a, AJ,   #
#                               BV, CDR3b, BJ)                       #
#                  -a --append (Append chains with constant regions  #
#                  -r --organism (Updated Organism - need for append #
#                                 Default "homo sapiens")            #
#                  -y --alpha (Alpha chain: constant region altern.  #
#                              Default "TRAC*01")                    #
#                  -z --beta (Beta chain: constant region altern.    #
#                             Default "TRBC1*01")                    #
#                  -m --omission COMMA SEP LIST (Characters to avoid #
#                                in header id for fasta file)        #
#                  -v --verbose (Show alignments being produced)     #
#                  -l --silent (Mute even poor alignments)           #
######################################################################


import argparse
import pandas
from Bio import Align


#################
#     Global    #
#################
gene_dic = {}
verbose = False
silent = False


#################
#    Methods    #
#################
def create_gene_dic(fasta_file, organism):
    """
    Populate the gene dictionary from IMGT data.

    Parameters
    __________
    fasta_file : str
        Fasta file containing IMGT data of gene segments of TCRs
    organism : str
        Organism name to pull from IMGT data file

    Returns
    _______
    gene_dic : dict
        Gene segment name -- AA sequence
    """
    with open(fasta_file, "r") as family:
        temp_fam = ""  # Holds current gene segment name
        for line in family:
            if line[0] == ">" and line.split("|")[2].upper() == organism.upper():  # Find header
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


def get_tcr_info(file, sheet, positions):
    """
    Create dictionary of components of TCR chains from table

    Parameters
    __________
    file : str
        .xlsx or .csv file containing AV, AJ, BV, BJ, CDR3A, CDR3B
    sheet : str
        optional parameter for .xlsx file to collect data from a secondary sheet (default="Sheet1")
    positions : str
        3: cloneID, 4: AV, 5: CDR3A, 7: AJ, 8: BV, 9: CDR3B, 11: BJ

    Returns
    _______
    tcr_dic : dict
        clone ID Positions: 0 - TRAV, 1 - CDR3A, 2 - TRAJ, 3 - TRBV, 4 - CDR3B, 5 - TRBJ, clone ID - 6
    """
    tcr_dic = {}  # Creates output dictionary
    clone_id_count = {}  # Keeps track of repeat clone ids
    previous_id = ""
    # Open if xlsx
    if file.endswith(".xlsx"):
        if sheet == "0":  # If no sheet name was provided, set default to primary sheet
            sheet = 0
        df = pandas.read_excel(file, sheet_name=sheet, engine='openpyxl')
    # Open if csv
    elif file.endswith(".csv"):
        df = pandas.read_csv(file)
    # if not clone_id provided
    if len(positions) == 6:
        positions.insert(0, -1)
        df["clone_ID"] = ""  # Add column to end of dataframe for the new clone_ids, referenced in position -1
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
                if not silent:
                    print(f"Skipped due to gene segment not found:\n{value}")
    return tcr_dic


def confirm_segment(seg_list):
    """
    Loop through gene segments of a clone id and return if a segment is not found in gene_dic

    Parameters
    __________
    seg_list : lst
        List of gene segments w/ clone id position 0

    Returns
    _______
    confirm segment : boolean
        True if gene segment discovered in gene_dic
    """
    for seg in seg_list:
        if seg not in gene_dic.keys():
            return False
    return True


def make_tcr_seq(tcr_parts):
    """
    Sends tcr parts into align_overlap to create resulting Alpha and Beta chain sequences.

    Parameters
    __________
    tcr_parts : dict
        dictionary containing each aa segment of TCR as constructed from get_tcr_info

    Returns
    _______
    output_segs : dict
        a dictionary with clone ID key and a list of Alpha and Beta chain sequences.
    """
    output_seqs = {}
    # For each clone-id
    for each in tcr_parts:
        # Alignment construction for TRAV, and CDR3A
        if verbose and not silent:
            print("\n\n\n>" + each)
            print("ALPHA: ")
            print("V: \n" + gene_dic[tcr_parts[each][0]])
            print("CDR: \n" + tcr_parts[each][1])
        temp_seq_alpha = align_overlap(gene_dic[tcr_parts[each][0]], tcr_parts[each][1], "V", each)
        if verbose and not silent:
            print("V-CDR: \n" + temp_seq_alpha)
        # Alignment construction for TRAV + CDR3A and TRAJ
        temp_seq_alpha = align_overlap(temp_seq_alpha, gene_dic[tcr_parts[each][2]], "J", each)
        if verbose and not silent:
            print("Full: \n" + temp_seq_alpha)
        # # Alignment construction for TRBV, and CDR3B
        if verbose and not silent:
            print("\nBETA: ")
            print("V: \n" + gene_dic[tcr_parts[each][3]])
            print("CDR: \n" + tcr_parts[each][4])
        temp_seq_beta = align_overlap(gene_dic[tcr_parts[each][3]], tcr_parts[each][4], "V", each)
        if verbose and not silent:
            print("V-CDR: \n" + temp_seq_beta)
        # Alignment construction for TRBV + CDR3B and TRBJ
        temp_seq_beta = align_overlap(temp_seq_beta, gene_dic[tcr_parts[each][5]], "J", each)
        if verbose and not silent:
            print("Full: \n" + temp_seq_beta)
        parts = tcr_parts[each]
        # 0: Alpha_chain, 1: Beta chain, 2: aV segment, 3: CDR3a, 4: aJ segment, 5: bV segment, 6: CDR3b, 7: bJ segment
        # 8: aV seq, 9: aJ seq, 10: bV seq, 11: bJ seq
        output_seqs[each] = [temp_seq_alpha, temp_seq_beta, parts[0], parts[1], parts[2], parts[3], parts[4], parts[5],
                             gene_dic[tcr_parts[each][0]], gene_dic[tcr_parts[each][2]], gene_dic[tcr_parts[each][3]],
                             gene_dic[tcr_parts[each][5]]]
    return output_seqs


# Method: 05 - align_overlap()
# Input:
#   front - start of seq
#   end - end of seq
#   part - If overlapping V or J segment
# Goal: Looks for overlap between front and end AAs besides when single AA overlap
# Output: Resulting overlapped segments -- if == J: full TCR sequence
def align_overlap(front, end, part, clone_id):
    """
    Looks for overlap between front and end AAs besides when single AA overlap

    Parameters
    __________
    front : str
        start of seq
    end : str
        end of seq
    part : str
        "V" or "J" depending on what segments are being overlapped to the cdr
    clone_id : str
        reference to what constant region to append to end of variable region sequence

    Returns
    _______
    output : str
        Resulting overlapped segments -- if == J and clone_id: full TCR sequence
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    best_align = []
    # Finding shortest alignment with best score to append regions
    # overlap cut modified for v or j based on testing
    if part == "V":
        # Alignment Weights
        aligner.match_score = 2.5
        aligner.mismatch_score = -1
        aligner.gap_score = -1.5
        aligner.extend_gap_score = -0.5
        # Perform alignment
        # When an alignment conforms to standard of at least one of the last 3 aa's of the V segment included in align
        good = False
        if len(end) > 7:
            overlap = (len(end) * -1) + 5  # Will decrease overlap if amino acids don't align properly (add from len)
        else:
            overlap = (len(end) * -1)  # If short CDR provided
        while not good:
            temp_align = aligner.align(front[overlap:], end)
            # Reference first alignment
            align_ = temp_align[0]
            # Collect alignment info
            front_seq = temp_align.alignment.format().split("\n")[0]
            pairs_info = temp_align.alignment.format().split("\n")[1]
            end_seq = temp_align.alignment.format().split("\n")[2]
            if verbose and not silent:
                print("\n" + front_seq)
                print(pairs_info)
                print(end_seq + "\n")
            # If last match is within the range needed
            if end_seq[0] == " " and front_seq[0] != " ":
                best_align = []
                best_align.append(overlap)  # overlap
                best_align.append(temp_align.alignment.format().split("\n")[1].count(" "))  # Count start position
                good = True
            else:
                if pairs_info.count("|") < 2:  # Only a single aa match - force append CDR3 and warn
                    best_align = "BAD"
                    good = True
                    if not silent:
                        print("Poor alignment possible: " + clone_id)
                elif len(front_seq) != 2:  # Bump overlap (trim front_seq) by one if not a good alignment
                    overlap += 1
                    if verbose and not silent:
                        full_terminal_text = "========================================================================"\
                                             "======================================================="
                        print(full_terminal_text)
    elif part == "J":
        offset = 6  # Helps with alignments to avoid too spaced out alignment
        # Alignment Weights
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.gap_score = -1
        aligner.extend_gap_score = -0.5
        # Perform alignment
        temp_align = aligner.align(front[(len(end) - offset) * -1:], end)
        align_ = temp_align[0]  # Reference first align
        # Collect alignment info to print
        front_seq = temp_align.alignment.format().split("\n")[0]
        pairs_info = temp_align.alignment.format().split("\n")[1]
        end_seq = temp_align.alignment.format().split("\n")[2]
        if verbose and not silent:
            print("\n" + front_seq)
            print(pairs_info)
            print(end_seq + "\n")
        # Collected needed information from alignment
        start_seq = temp_align.alignment.format().split("\n")[0]
        end_seq = temp_align.alignment.format().split("\n")[2]
        end_pos = len(temp_align.alignment.format().split("\n")[1])
        best_align = [start_seq, end_seq, end_pos]
    if part == "V":
        if best_align != "BAD":
            end_front = len(front) + best_align[0] + best_align[1]
            output = front[:end_front] + end
        else:
            output = front + end
        return output
    elif part == "J":  # If joining region overlap
        output = front + best_align[1][best_align[2]:]  # Overlaps front with remainder of best_align end seq
        return output


def if_gap(front, end, region, gap_count):  # Check for cdr len in previous method, send in if V or J region.
    """
    Determine if there are gaps in the alignments and where in the sequence it's located

    Parameters
    __________
    front : str
        start of seq
    end : str
        end of seq
    region : str
        "V" or "J" depending on what segments are being overlapped to the cdr
    gap_count : int
        location of gap

    Returns
    _______
    seq : str
        Adjusted aligned sequence
    """
    temp_seq = ""
    if region == "V":  # When the variable region and cdr are being overlapped.
        for letter in front:
            if letter != "-":
                temp_seq += letter
            else:
                return temp_seq + end[len(temp_seq):]
    elif region == "J":  # When the joining region and cdr are being overlapped.
        return front + end[gap_count:]


def append_constant(tcr_seq_dic, organism, alpha, beta):
    """
    Take in the variable region sequence and append the constant region based on the organism name provided

    Parameters
    __________
    tcr_seq_dic : dict
        dictionary of tcr components clone_id: alpha, beta, [other parts]
    organism : str
        organism name provided by user
    alpha : str
        alpha chain constant region id
    beta : str
        beta chain constant region id

    Return
    ______
    tcr_seq_dic : dict
        updated dictionary with all tcr sequences appended with constant regions
    """
    constant_dic = create_constant_dic(organism)
    for tcr_id in tcr_seq_dic:
        for segment in constant_dic:
            if constant_dic[segment][0].startswith(alpha) and constant_dic[segment][1] == "EX1":
                tcr_seq_dic[tcr_id][0] = tcr_seq_dic[tcr_id][0] + constant_dic[segment][2].replace("X", "")
            if constant_dic[segment][0].startswith(beta) and constant_dic[segment][1] == "EX1":
                tcr_seq_dic[tcr_id][1] = tcr_seq_dic[tcr_id][1] + constant_dic[segment][2].replace("X", "")
    return tcr_seq_dic


def create_constant_dic(organism):
    """
    Populate the gene dictionary from IMGT data.

    Parameters
    __________
    organism : str
        name of organism provided by user

    Returns
    _______
    constant_dic : dict
        dictionary of the constant regions for the organism selected
    """
    constant_dic = {}
    with open("constant_regions.fasta", "r") as f:
        temp_imgt_id = ""
        flag = False  # Flag to know when to collect sequence information for dic
        for line in f:
            if line[0] == ">":  # Determine if header
                if line.split("|")[2].upper() == organism.upper():  # Determine if correct organism
                    header = line.split("|")
                    temp_imgt_id = header[0][1:] + "-" + header[4]
                    segment_name = header[1]  # ex. TRAC*01
                    exon_num = header[4]  # ex. EX1
                    constant_dic[temp_imgt_id] = [segment_name, exon_num]  # prepares location in dictionary
                    flag = True
                else:
                    flag = False
            elif flag:
                # Collecting sequence data from fasta file per header
                if len(constant_dic[temp_imgt_id]) == 2:  # Catches when multiple lines of seq
                    if line.__contains__('*'):
                        constant_dic[temp_imgt_id].append(line[:-1].replace('*', ""))
                    elif line.__contains__('X'):
                        constant_dic[temp_imgt_id].append(line[:-1].replace('X', ""))
                    else:
                        constant_dic[temp_imgt_id].append(line[:-1])
                elif len(constant_dic[temp_imgt_id]) == 3:
                    if line.__contains__('*'):
                        constant_dic[temp_imgt_id][2] += (line[:-1].replace('*', ""))
                    elif line.__contains__('X'):
                        constant_dic[temp_imgt_id][2] += (line[:-1].replace('X', ""))
                    else:
                        constant_dic[temp_imgt_id][2] += (line[:-1])
    return constant_dic


def make_info_table(tcr_dic, file_in="translated_seq.csv"):
    """
    Create information table to return to the use with details about the constructed TCR sequences

    Parameters
    __________
    tcr_dic : dict
        resulting clone ids and tcr information
    file_in : str
        optional file name for information table
    """
    if "/" in file_in:
        file_in = file_in.split("/")[-1]
    with open(file_in, "w") as f:  # Open file, name is set by default but can be changed
        # Header
        f.write("clone ID,TRAV,CDR3 Alpha,TRAJ,TRBV,CDR3 Beta,TRBJ,TRAV_seq,TRAJ_seq,TRBV_seq,TRBJ_seq,ALPHA,BETA\n")
        for each in tcr_dic:  # Loop through each clone ID
            info = tcr_dic[each][2:]
            chains = tcr_dic[each][:2]
            f.write(each + "," + ",".join(info) + "," + ",".join(chains) + "\n")


def make_fasta_files(alpha_file, beta_file, tcr_dic, character_list):
    """
    Create alpha.fasta and beta.fasta - containing alpha and beta chain sequences with id as header

    Parameters
    __________
    alpha_file : str
        name of alpha chain file
    beta_file : str
        name of beta chain file
    tcr_dic : dict
        dictionary of all tcrs
    character_list : str
        list of characters to not include in fasta header id
    """
    with open(alpha_file, "w") as file_a:
        with open(beta_file, "w") as file_b:
            for clone_id in tcr_dic:
                clone_id_print = clone_id
                # Remove characters not allowed in clone_id
                if character_list != "":
                    for char in character_list.split(","):
                        if char in clone_id:
                            clone_id_print = clone_id.replace(char, "")
                tcr_alpha_chain = tcr_dic[clone_id][0]
                tcr_beta_chain = tcr_dic[clone_id][1]
                count_1 = 1
                count_2 = 1
                if tcr_alpha_chain != '':
                    file_a.write('>' + clone_id_print + "\n")
                    for aa in tcr_alpha_chain:
                        if count_1 % 81 != 0:
                            file_a.write(aa)
                            count_1 += 1
                        else:
                            file_a.write('\n' + aa)
                            count_1 = 2
                    file_a.write('\n')
                if tcr_beta_chain != '':
                    file_b.write('>' + clone_id_print + '\n')
                    for aa in tcr_beta_chain:
                        if count_2 % 81 != 0:
                            file_b.write(aa)
                            count_2 += 1
                        else:
                            file_b.write('\n' + aa)
                            count_2 = 2
                    file_b.write('\n')


####################
#     Controls     #
####################
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("single_cell_table", help="Xlsx or Csv of single cell t-cell data", type=str)
    parser.add_argument("-s", "--sheet", help="Sheet name in table", default="0", type=str)
    parser.add_argument("-f", "--fasta", help="(Output) Fasta file format of results", default=False,
                        action="store_true")
    parser.add_argument("-t", "--information", help="(Output) Creates information table of constructed TCR seq."
                        , default=True, action="store_true")
    parser.add_argument("-g", "--genefamily",
                        help='Replacement gene family segment file (if not homo sapiens) | Default = "family_seq.fasta"',
                        default="family_seq.fasta", type=str)
    parser.add_argument("-c", "--columns",
                        help="Positions of gene segments: Clone ID,AV,CDR3a,AJ,BV,CDR3b,BJ | Can exclude Clone ID",
                        default="0,1,2,3,4,5,6", type=str)  # Our data 3,4,5,7,8,9,11
    parser.add_argument("-a", "--append", help="Append each chain with constant region", default=False,
                        action="store_true")
    parser.add_argument("-r", "--organism", help="Updated Organism",
                        default="homo sapiens", type=str)
    parser.add_argument("-y", "--alpha", help='Alpha Chain: Constant region alternative | Deafult = "TRAC*01"',
                        type=str, default="TRAC*01")
    parser.add_argument("-z", "--beta", help='Beta Chain: Constant region alternative | Deafult = "TRBC1*01"', type=str,
                        default="TRBC1*01")
    parser.add_argument("-m", "--omission", help="Characters to avoid in header id for fasta file submission | comma"\
                        "sep.", type=str, default="")
    parser.add_argument("-v", "--verbose", help="Show alignments being produced", default=False, action="store_true")
    parser.add_argument("-l", "--silent", help="Mute even poor alignments", default=False, action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()  # collect arguments
    if args.verbose:
        global verbose
        verbose = True
    if args.silent:
        global silent
        silent = True
    # create dictionary for program to pull gene family information from
    create_gene_dic(args.genefamily, args.organism, args.alpha, args.beta)
    # create tcr dictionary containing each segment
    tcr_dic = get_tcr_info(args.single_cell_table, args.sheet, [int(i) for i in args.columns.split(",")])
    # create full tcr sequences
    tcr_seq_dic = make_tcr_seq(tcr_dic)
    if args.append:
        tcr_seq_dic = append_constant(tcr_seq_dic, args.organism)
    # Create information table with full TCR sequences
    if args.fasta:
        make_fasta_files("alpha.fasta", "beta.fasta", tcr_seq_dic, args.omission)
    elif args.information:
        make_info_table(tcr_seq_dic, ".".join(args.single_cell_table.split(".")[:-1]) + "_Results.csv")


if __name__ == '__main__':
    main()
