#!/usr/bin/python3

######################################################################
# PostCoupler.py -- A component of TRain                             #
# Copyright: Austin Seamann & Dario Ghersi                           #
# Version: 0.1                                                       #
# Last Updated: January 12th, 2022                                   #
# Goal: Match chain labeling and numbering of a reference PDB file   #
#                                                                    #
# Positional argument: PDB (PDB file being renumbered                #
# Named arguments: -r --reference (Reference PDB file being used to  #
#                                 match numbering)                   #
#                  -n --name (Customized name for output file)       #
#                  -c --chains (Chains that will be renumbered from  #
#                              reference)                            #
#                  -u --custom (Custom numbering for smalled chain   #
#                              | e.g. 9,11,2,4,5 | Has to match num  #
#                              of aa to numbers provided)            #
######################################################################


import argparse
from Bio import pairwise2


class MatchNumber:
    #################
    #    Methods    #
    #################
    def match_number(self, custom):
        """
        Align chains and for amino acids that overlap, renumber them to the reference chains
        numbering. Update them to PDB.

        Parameters
        __________
        custom : str
            comma sep. str of numbering for shorest chain if provided by user
        """
        # Constants
        POSSEQ = [22, 26]
        CHAINID = 21
        AA = [17, 20]

        ref_chains = self.get_chains(self.reference)
        ref_seqs = {}  # Dictionary of seq of ref_seqs chain id: seq
        target_chains = self.get_chains(self.target)
        target_seqs = {}  # Dictionary containing chain id: seq
        ref_aa = {}  # Dictionary containing aa numbering of ref.

        for chain in target_chains:
            ref_seqs[chain] = self.get_amino_acid_on_chain(self.reference, chain)
            target_seqs[chain] = self.get_amino_acid_on_chain(self.target, chain)
            if chain in ref_chains:
                ref_aa[chain] = {}  # Dictionary for num aa: aa of target.
        with open(self.reference, "r") as r:
            for line in r:
                if line[0:6] == 'ATOM  ':  # Only references atoms
                    if line[CHAINID] in ref_aa.keys():
                        ref_aa[line[CHAINID]][int(line[POSSEQ[0]:POSSEQ[1]])] = line[AA[0]:AA[1]]
        aligns = {}  # alignments to determine start of chain number if they don't start at the same point
        for chain in target_chains:
            aligns[chain] = pairwise2.align.globalms(target_seqs[chain], ref_seqs[chain], 2, -1, -2, -.5,
                                                     penalize_end_gaps=(False, False), one_alignment_only=True)
            # If there is a gap at the start of the target seq, removes those positions from ref_aa until it matches
            start_of_chain = self.find_gap(aligns[chain][0])
            while start_of_chain != 0:
                ref_aa[chain].pop(ref_aa[chain].key()[0])
                start_of_chain -= 1
        shortest_num = 10000
        shortest_chain = ""
        # For custom numbering, determining shortest chain
        for chain in target_seqs:
            if len(target_seqs[chain]) < shortest_num:
                shortest_chain = chain
                shortest_num = len(target_seqs[chain])
        if custom != "...":
            custom_num = custom.split(",")
            ref_aa[shortest_chain] = {}
            for aa in target_seqs[shortest_chain]:
                ref_aa[shortest_chain][custom_num.pop(0)] = self.one_to_three(aa)
        with open(self.target, "r") as t:
            with open(self.name, "w") as w:
                atom_count = 1  # No matter what, starts atom count at 1
                aa_num = -10000
                current_chain = ""
                for line in t:
                    if line[0:6] == 'ATOM  ':
                        if aa_num < 0:  # Checks to see if we're just starting
                            aa_num = int(line[POSSEQ[0]: POSSEQ[1]])  # Tell me current amino acid number
                            current_chain = line[CHAINID]  # Tells us the starting chain
                        if current_chain != line[CHAINID]:  # Checks to see if we're on a new chain
                            ref_aa.pop(current_chain)  # Removes previous chain from list
                            current_chain = line[CHAINID]  # Updates current chain
                            aa_num = int(line[POSSEQ[0]: POSSEQ[1]])  # Updates current aa count
                        if aa_num != int(line[POSSEQ[0]: POSSEQ[1]]):  # Tells us if we're on a new aa
                            ref_aa[line[CHAINID]].pop(list(ref_aa[line[CHAINID]])[0])  # Removes last counted AA
                            aa_num = int(line[POSSEQ[0]: POSSEQ[1]])
                        temp_line = ""
                        temp_line += line[0:6]  # Adds header
                        temp_line += str(atom_count).rjust(5)  # Adds atom count
                        atom_count += 1
                        if self.chains.__contains__(line[CHAINID]):  # Skips chains that we're not renumbering
                            temp_line += line[11:POSSEQ[0]]  # Replaces gap from num to AA num
                            next_num = list(ref_aa[line[CHAINID]])[0]
                            temp_line += str(next_num).rjust(4)
                            temp_line += line[POSSEQ[1]:]
                            w.write(temp_line)
                        else:  # If we're skipping a chain we just write it but we keep atom numbering
                            temp_line += line[11:]
                            w.write(temp_line)
                    else:
                        w.write(line)

    def find_gap(self, seq):
        """
        Counts how may gaps until start of seq and reports position

        Parameters
        __________
        seq : str
            sequence string to be search for gap character ("-")

        Returns
        _______
        count : int
            position in which first gap character appears
        """
        count = 0
        for letter in seq:
            if letter == "-":
                count += 1
            else:
                break
        return count

    def get_chains(self, pdb):
        """
        Collect chains found in PDB file

        Parameters
        __________
        pdb : str
            location of pdb file

        Returns
        _______
        chains : list
            list of all chains in PDB file
        """
        chains = []
        with open(pdb, 'r') as file:
            for line in file:
                if line[0:6] == 'ATOM  ':
                    if not chains.__contains__(line[21]):
                        chains.append(line[21])
        return chains

    # Returns a string of amino acids in a specific chain as a string in single letter notation
    # INPUT: PDB name, ChainID
    def get_amino_acid_on_chain(self, pdb, chain):
        output = ''
        count = 0
        flag = True
        with open(pdb, 'r') as file:
            for line in file:
                if line[0:6] == 'ATOM  ':
                    if line[21] == chain:
                        if flag:
                            count = int(line[23:26])
                            flag = False
                        if count == int(line[23:26]):
                            if line[16] != 'B':
                                output += self.three_to_one(line[17:20])
                                count += 1
                        elif count < int(line[23:26]):
                            count = int(line[23:26])
        return output

    # Converts three letter AA to single letter abbreviation
    def three_to_one(self, three):
        translate = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'ASX': 'B', 'CYS': 'C', 'GLU': 'E',
            'GLN': 'Q', 'GLX': 'Z', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
            'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y',
            'VAL': 'V'
        }
        for key in translate:
            if three.upper() == key:
                return translate[key]

    # Converts single letter abbreviation to three letter AA
    def one_to_three(self, one):
        translate = {
            'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'B': 'ASX', 'C': 'CYS', 'E': 'GLUE',
            'Q': 'GLN', 'Z': 'GLX', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
            'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR',
            'V': 'VAL'
        }
        for key in translate:
            if one.upper() == key:
                return translate[key]

    # Set output pdb name
    # INPUT: name
    def set_name(self, inname):
        self.name = inname

    def set_chains(self, chains_in):
        self.chains = chains_in.upper()

    # Initializes class with target and reference
    def __init__(self, target, reference):
        if target[-4:] == ".pdb":
            self.target = target
        if reference[-4:] == ".pdb":
            self.reference = reference
        self.chains = self.get_chains(self.target)
        self.name = target[:-4] + "_orinum.pdb"


####################
#     Controls     #
####################
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="PDB file being renumbered", type=str)
    parser.add_argument("-r", "--reference",
                        help="Reference PDB file being used to match numbering",
                        type=str)
    parser.add_argument("-n", "--name",
                        help="Customized name for output file", type=str)
    parser.add_argument("-c", "--chains",
                        help="Chains that will be renumbered from ref.",
                        type=str)
    parser.add_argument("-u", "--custom", help="Custom numbering for smallest chain | e.g. 9,11,2,4,5 | Has to match "\
                        "num of aa to numbers provided", type=str,
                        default="...")
    return parser.parse_args()


def main():
    args = parse_args()
    match = MatchNumber(args.pdb, args.reference)
    if args.name:
        match.set_name(args.name)
    if args.chains:
        match.set_chains(args.chains)
    match.match_number(args.custom)


if __name__ == '__main__':
    main()
