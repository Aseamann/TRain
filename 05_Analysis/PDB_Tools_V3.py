#!/usr/bin/python3

######################################################################
# PDB_Tools_V3.py -- A component of TRain                            #
# Copyright: Austin Seamann, Dario Ghersi, and Ryan Ehrlich          #
# Version: 0.1                                                       #
# Last Updated: June 26th, 2022                                      #
# Goal: Be a tool that is resourceful for working with PDB files in  #
#       general and for the specific use of T-cell Receptor          #
#       PDB files.                                                   #
#                                                                    #
# Positional argument: PDB File in Use                               #
# Named arguments: SEE -h LIST IN CONSOLE OR VISIT BOTTOM OF FILE    #
######################################################################


import argparse
from math import sqrt, pow
import statistics
import numpy as np
from sklearn.decomposition import PCA
from scipy.spatial.transform import Rotation
from math import sqrt
from Bio import Align
from Bio.Align import substitution_matrices
import Bio.PDB
import os

#################
#     Global    #
#################
# PDB indices
AA = [19, 70]
POSX = [30, 38]
POSY = [38, 46]
POSZ = [46, 54]
POSRES = [17, 20]
POSSEQ = [22, 26]
POSCHAIN = 21
ChainID = 11


#################
#    Methods    #
#################
class PdbTools3:
    # initialize PdbTools
    def __init__(self, file="..."):
        """
        Initialize PdbTools

        Parameters
        __________
        file : str
            PDB file in use
        """
        self.file_name = file
        self.test_list = {}

    def set_file_name(self, file_name_in):
        """
        Method for changing PDB file in use

        Parameters
        __________
        file_name_in : str
            name of new PDB to be in use

        """
        self.file_name = file_name_in

    def get_file_name(self):
        """
        Returns file name currently in use

        Returns
        _______
        file_name : str
            name of PDB file
        """
        return self.file_name

    def get_pdb_id(self):
        """
        Returns the PDB ID of file

        Returns
        _______
        PDB id based on the text in the header of PDB
        """
        with open(self.file_name, 'r') as file:
            output = file.readline()
            return output[62:66].lower()

    # Returns a list of all chains in PDB file
    def get_chains(self):
        """
        Returns a list of all chains in PDB file

        Returns
        _______
        chains : list
            List of chains contained in PDB file
        """
        chains = []
        with open(self.file_name, 'r') as file:
            for line in file:
                if line[0:6] == 'ATOM  ':
                    if not chains.__contains__(line[21]):
                        chains.append(line[21])
        return chains

    def get_resolution(self):
        """
        Reads the PDB for the resolution of the file

        Returns
        -------
        output : float
            Resolution of file if contained in PDB file
        """
        value = ''
        with open(self.file_name, 'r') as file:
            flag = False
            for line in file:
                if line[0:6] == 'REMARK' and not flag:
                    temp = line[6:10]
                    if temp == '   2':
                        flag = True
                elif line[0:6] == 'REMARK' and flag:
                    value += line[26:29]
                    break
        output = float(value)
        return output

    def get_amino_acid_on_chain(self, chain):
        """
        Returns a string of amino acids in a specific chain as a string in single letter notation

        Parameters
        ----------
        chain : str

        Returns
        -------
        output : str
            String of aa's in single letter formatting
        """
        output = ''
        count = 0
        flag = True
        with open(self.file_name, 'r') as file:
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

    def first_atom_on_chain(self, chain):
        """
        Returns a dictionary for the first atom of a chain.

        Parameters
        ----------
        chain : str

        Returns
        -------
        atom : dict
            Contains elements of the first atom in a chain
        """
        with open(self.file_name, 'r') as file:
            for line in file:
                if line[0:6] == 'ATOM  ':
                    if line[21] == chain.upper() and len(line) >= 76:
                        atom = {'atom_num': int(line[6:11]), 'atom_id': line[13:16].strip(),
                                'atom_comp_id': line[17:20],
                                'chain_id': line[21], 'comp_num': int(line[22:26]), 'X': float(line[31:38]),
                                'Y': float(line[38:46]), 'Z': float(line[46:54]), 'occupancy': float(line[55:60]),
                                'B_iso_or_equiv': float(line[60:66]), 'atom_type': line[77]}
                        return atom
                    elif line[21] == chain.upper() and len(line) <= 76:
                        atom = {'atom_num': int(line[6:11]), 'atom_id': line[13:16].strip(),
                                'atom_comp_id': line[17:20],
                                'chain_id': line[21], 'comp_num': int(line[22:26]), 'X': float(line[31:38]),
                                'Y': float(line[38:46]), 'Z': float(line[46:54]), 'occupancy': float(line[55:60])}
                        return atom

    def three_to_one(self, three):
        """
        Converts three letter AA to single letter abbreviation

        Parameters
        ----------
        three : str
            Three letter abbreviation of amino acid

        Returns
        -------
        Converted amino acid format
        """
        translate = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'ASX': 'B', 'CYS': 'C', 'GLU': 'E',
            'GLN': 'Q', 'GLX': 'Z', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
            'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y',
            'VAL': 'V'
        }
        for key in translate:
            if three.upper() == key:
                return translate[key]

    def get_atom(self, atom_num):
        """
        Returns a dictionary with elements related to specific atom number entered

        Parameters
        ----------
        atom_num : int
            Integer value for atom in file

        Returns
        -------
        atom : dict
            Dictionary containing information for atom from PDB file
        """
        with open(self.file_name, 'r') as file:
            for line in file:
                if line[0:6] == 'ATOM  ':
                    if int(line[6:11]) == atom_num and len(line) >= 76:
                        atom = {'atom_num': int(line[6:11]), 'atom_id': line[13:16].strip(),
                                'atom_comp_id': line[17:20],
                                'chain_id': line[21], 'comp_num': int(line[22:26]), 'X': float(line[31:38]),
                                'Y': float(line[38:46]), 'Z': float(line[46:54]), 'occupancy': float(line[55:60]),
                                'B_iso_or_equiv': float(line[60:66]), 'atom_type': line[77]}
                        return atom
                    elif int(line[6:11]) == atom_num and len(line) <= 76:
                        atom = {'atom_num': int(line[6:11]), 'atom_id': line[13:16].strip(),
                                'atom_comp_id': line[17:20],
                                'chain_id': line[21], 'comp_num': int(line[22:26]), 'X': float(line[31:38]),
                                'Y': float(line[38:46]), 'Z': float(line[46:54]), 'occupancy': float(line[55:60])}
                        return atom

    def euclidean_of_atoms(self, atom_num_1, atom_num_2):
        """
        Returns the Euclidean distance between two atoms based with atom_id being sent in as the parameter

        Parameters
        ----------
        atom_num_1 : int
        atom_num_2 : int

        Returns
        -------
        euclidean_distance : float
            Distance between two atoms in angstroms
        """
        atom_1 = self.get_atom(atom_num_1)
        atom_2 = self.get_atom(atom_num_2)
        euclidean_distance = sqrt((atom_2['X'] - atom_1['X'])**2 + (atom_2['Y'] - atom_1['Y'])**2
                                  + (atom_2['Z'] - atom_1['Z'])**2)
        return euclidean_distance

    # Collect the atoms from an inputted chain. Provides all values in PDB file
    def get_atoms_on_chain(self, chain):
        """
        Collect the atoms from an inputted chain. Provides all values in PDB file

        Parameters
        ----------
        chain : str

        Returns
        -------
        atoms : list
            List of atoms based on the chain submitted
        """
        atoms = []
        with open(self.file_name, 'r') as file:
            for line in file:
                if line[0:6] == 'ATOM  ':
                    if line[21] == chain and len(line) >= 76:
                        atoms.append({'atom_num': int(line[6:11]), 'atom_id': line[13:16].strip(),
                                'atom_comp_id': line[17:20],
                                'chain_id': line[21], 'comp_num': int(line[22:26]), 'X': float(line[31:38]),
                                'Y': float(line[38:46]), 'Z': float(line[46:54]), 'occupancy': float(line[55:60]),
                                'B_iso_or_equiv': float(line[60:66]), 'atom_type': line[77]})
                    elif line[21] == chain and len(line) >= 76:
                        atoms.append({'atom_num': int(line[6:11]), 'atom_id': line[13:16].strip(),
                                'atom_comp_id': line[17:20],
                                'chain_id': line[21], 'comp_num': int(line[22:26]), 'X': float(line[31:38]),
                                'Y': float(line[38:46]), 'Z': float(line[46:54]), 'occupancy': float(line[55:60])})
        return atoms

    def rebuild_atom_line(self, atoms):
        """
        Rebuild the string of an atom based on atom dictionary element of PDB_Tools

        Parameters
        ----------
        atoms : int

        Returns
        -------
        output : str
        """
        output = ""
        for atom in atoms:
            # Reformat ATOM lines - Write this in separate method
            line = 'ATOM  ' + str(atom['atom_num']).rjust(5) + "  " + atom['atom_id'].ljust(3) + " "
            line += atom['atom_comp_id'] + " " + atom['chain_id'] + str(atom['comp_num']).rjust(4) + "    "
            # XYZ
            line += str(format(atom['X'], ".3f")).rjust(8) + str(format(atom['Y'], ".3f")).rjust(8)
            line += str(format(atom['Z'], ".3f")).rjust(8)
            # Post XYZ
            line += str(format(atom['occupancy'], ".2f")).rjust(6)
            line += str(format(atom['B_iso_or_equiv'], ".2f")).rjust(6) + "           "
            line += atom['atom_type'] + "\n"
            output += line
        return output

    def get_tcr_chains(self):
        """
        Returns alpha and beta chain IDs based on seq. alignment to 1a07 PDB entry chains. Confirms that it is a
        partnering TCR chain

        Returns
        -------
        result : str
        """
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
        result = {}
        # Hard coded peptide chains for alpha and beta elements of the TCR_file
        alpha_chain = [
            'KEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFTYRQYSGKSPELIMSIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSDSATYLCAVTTDSTGKLQFGAGT'\
            'QVVVTPDIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVATSNKSDFACANAFNNSIIPEDTFFPSPESS',
            'QKVTQTQTSISVMEKTTVTMDCVYETQDSSYFLFTYKQTASGEIVFLIRQDSYKKENATVGHYSLNFQKPKSSIGLIITATQIEDSAVYFCAMRGDYGGSGNKL'\
            'IFGTGTLLSVKP']
        beta_chain = [
            'NAGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSTYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQ'\
            'YFGPGTRLTVTEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSTTVNGKEVHSGVSTDPQPLKEQPALNDSRYALSSRLRVSATFTQNPRNHF'\
            'RCQVQFYGLSENDETTQDRAKPVTQIVSAEATGRAD',
            'VTLLEQNPRTRLVPRGQAVNLRCILKNSQYPTMSTYQQDLQKQLQTLFTLRSPGDKEVKSLPGADYLATRVTDTELRLQVANMSQGRTLYCTCSADRVGNTLYFG'\
            'EGSRLIV']
        chains = self.get_chains()
        tmp_alpha = []
        tmp_beta = []
        # Assume that alpha and beta chains are next to each other in PDB file
        for pos in range(0, len(chains)):
            score_alpha = aligner.align(self.get_amino_acid_on_chain(chains[pos]), alpha_chain[0])[0].score
            tmp_alpha.append([float(score_alpha), chains[pos]])
            # Position in front
            if pos + 1 in range(0, len(chains)):
                score_beta = aligner.align(self.get_amino_acid_on_chain(chains[pos + 1]), beta_chain[0])[0].score
                tmp_beta.append([float(score_beta), chains[pos + 1]])
            # Position behind
            if pos - 1 in range(0, len(chains)):
                score_beta = aligner.align(self.get_amino_acid_on_chain(chains[pos - 1]), beta_chain[0])[0].score
                tmp_beta.append([float(score_beta), chains[pos - 1]])
        alpha = sorted(tmp_alpha)
        beta = sorted(tmp_beta)
        result['ALPHA'] = alpha[-1][1]
        result['BETA'] = beta[-1][1]
        return result

    def get_tcr_amino_seq(self, tcr_type_in):
        """
        Returns the amino acid sequence of either 'ALPHA' or 'BETA' chain as single letter AA abbreviation

        Parameters
        ----------
        tcr_type_in : str
            'ALPHA' or 'BETA'

        Returns
        -------
        String of amino acid contained in desired chain.
        """
        tcr_dict = self.get_tcr_chains()
        for key in tcr_dict:
            if key == tcr_type_in:
                return self.get_amino_acid_on_chain(tcr_dict[key])

    def get_atoms_on_peptide(self):
        """
        Returns the atoms on the peptide, peptide is determined by smallest chain

        Returns
        -------
        atoms : list
        """
        atoms = []
        self.set_record_type('ATOM')
        file_atoms = self.record_report().splitlines()
        chain = self.get_peptide_chain()
        for line in file_atoms:
            if line[15] == chain:
                atoms.append(self.get_atom(int(line[0:5])))
        atoms.pop(len(atoms) - 1)
        return atoms

    def get_peptide_chain(self):
        """
        Returns the AA chain that is the peptide of the pMHC complex, based on the smallest chain in file

        Returns
        -------
        peptide : str
        """
        chains = self.get_chains()
        chain_dic = {}
        peptides = []
        for chain in chains:
            chain_dic[chain] = self.get_amino_acid_on_chain(chain)
        # Assumes peptide is shorter than 20 aa's
        for chain in chain_dic:
            if len(chain_dic[chain]) < 20:
                peptides.append(chain)
        # Peptide has to be within set distance of first AA in MHC
        mhc_1 = self.first_atom_on_chain(self.get_mhc_chain())
        position = 0
        while True:
            peptide_first = self.first_atom_on_chain(peptides[position])
            peptide_last = self.get_atoms_on_chain(peptides[position])[-1]
            # Distance between N-terminus of MHC and first AA in peptide
            distance1 = self.euclidean_of_atoms(mhc_1['atom_num'], peptide_first['atom_num'])
            # Distance between N-terminus of MHC and last AA in peptide
            distance2 = self.euclidean_of_atoms(mhc_1['atom_num'], peptide_last['atom_num'])
            # Determines if N-terminus is close enough to start or end of peptide being tested
            if distance1 <= 35 or distance2 <= 35:
                peptide = peptides[position]
                break
            else:
                position += 1
        return peptide

    def get_mhc_chain(self):
        """
        Returns the MHC chain based on a sequence alignment to a hardcoded MHC chain

        Returns
        -------
        mhc chain id
        """
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
        # Hard coded mhc chain
        mhc_chain = 'GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTV'\
                    'QRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMT'\
                    'HHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWE'
        chains = self.get_chains()
        tmp_mhc = []
        for chain in chains:
            score_mhc = aligner.align(self.get_amino_acid_on_chain(chain), mhc_chain)[0].score
            tmp_mhc.append([float(score_mhc), chain])
        mhc = sorted(tmp_mhc)
        return mhc[-1][1]

    def get_b2m_chain(self):
        """
        Return the Beta-2 Microglobulin - based on alignment to reference 1ao7 B2M

        Returns
        -------
        b2m chain id
        """
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
        # Hard coded b2m chain
        b2m_chain = 'MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYCTEFTPTEKDEYACRVNHVTLSQPCIVKWDRDM'
        chains = self.get_chains()
        tmp_b2m = []
        for chain in sorted(chains):
            score_b2m = aligner.align(self.get_amino_acid_on_chain(chain), b2m_chain)[0].score
            tmp_b2m.append([float(score_b2m), chain])
        b2m = sorted(tmp_b2m, reverse=True)
        high_score = b2m[0][0]  # highest align score
        b2ms = []
        for each in b2m:  # Tries to grab the first mhc in file
            if each[0] / high_score >= 0.98:  # similarity is above 98%
                b2ms.append(each[1])
        return sorted(b2ms)[0]

    def renumber_docking(self, rename="****"):
        """
        Return a TCR PDB with chains and atoms renumber for rosettadock.

        Parameters
        ----------
        rename : str
            Optional updated naming to PDB file being created
        """
        if rename != '****':  # Renames if given input, if not writes over file
            tcr = rename
        else:
            tcr = self.file_name
        atom_count = 1
        res_count = 1
        previous_res_count = 0
        previous_chain = ""
        flag_start_res = False
        file_save = ""
        with open(self.file_name, "r") as f:  # Reads in PDB
            for line in f:
                file_save += line
        with open(tcr, "w") as f1:  # Writes renumbered PDB
            for line in file_save.split("\n"):
                if line[0:6] == 'HEADER':
                    f1.write(line + "\n")
                if line[0:6] == 'ATOM  ':  # Only writes over atoms
                    if line[16] != 'B' and line[26] == ' ':  # Don't allow secondary atoms
                        if line[21] != previous_chain and previous_chain != "":  # Writes TER at end of chain
                            f1.write("TER\n")
                        if not flag_start_res:
                            previous_res_count = int(line[22:26])
                            flag_start_res = True
                        if int(line[22:26]) != previous_res_count:
                            previous_res_count = int(line[22:26])
                            res_count += 1
                        line = line[:22] + str(res_count).rjust(4) + line[26:]  # Replaces res count
                        line = line[:6] + str(atom_count).rjust(5) + line[11:]  # Replaces atom count
                        if line[16] == 'A':
                            line = line[:16] + ' ' + line[17:]
                        previous_chain = line[21]  # Updates previous chain for adding TER
                        atom_count += 1
                        f1.write(line + '\n')
            f1.write("END\n")

    def clean_tcr(self, dir_start='****'):
        """
        Returns a reformatted PDB with just the TCR atom cord. and has relabeled chains with ALPHA = A; BETA = B

        Parameters
        ----------
        dir_start : str
            Choose what directory to save PDB
        """
        if dir_start != '****':
            tcr = dir_start + '%s_tcr.pdb' % (self.get_pdb_id())
        else:
            tcr = '%s.pdb' % (self.get_pdb_id())
        tcr_list = self.get_tcr_chains()
        atom_count = 0
        flag = False
        output = []
        with open(self.file_name) as f:
            for line in f:
                if line[0:6] == 'HEADER':
                    output.append(line)
                    flag = True
                if flag:
                    output.append('EXPDTA    THEORETICAL MODEL    CLEAN TCR ALPHA:A BETA:B\n')
                    flag = False
                if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                    if line[16] != 'B':
                        if line[21] == tcr_list.get('ALPHA') or line[21] == tcr_list.get('BETA'):
                            if line[21] == tcr_list.get('ALPHA'):
                                line = line[:21] + 'A' + line[22:]
                            elif line[21] == tcr_list.get('BETA'):
                                line = line[:21] + 'B' + line[22:]
                            num = line[6:11]
                            atom_count += 1
                            if line[16] == 'A':
                                line = line[:16] + ' ' + line[17:]
                            output.append(line.replace(num, str(atom_count).rjust(5), 1))
        with open(tcr, 'w+') as f1:
            for line in output:
                f1.write(line)

    def clean_tcr_count_trim(self, dir_start='****'):
        """
        Returns a reformatted PDB with just the TCR atom cord. and has relabeled chains with ALPHA = A; BETA = B
        Also uses trimmed TCR seq.
        New count on amino acids so each chain starts at 1

        Parameters
        ----------
        dir_start : str
            Choose what directory to save PDB
        """
        if dir_start != '****':
            tcr = dir_start + '%s_tcr.pdb' % (self.get_pdb_id())
        else:
            tcr = self.get_pdb_id() + ".pdb"
        tcr_list = self.get_tcr_chains()
        atom_count = 0
        flag = False
        flag_a = False
        flag_b = False
        res_alpha_count = 1
        alpha_cut = 107
        res_beta_count = 1
        beta_cut = 113
        previous_count_a = -1
        previous_count_b = -1
        res_count = 0
        output = []
        with open(self.file_name) as f:
            for line in f:
                if line[0:6] == 'HEADER':
                    output.append(line)
                    flag = True
                if flag:
                    output.append('EXPDTA    THEORETICAL MODEL    CLEAN TCR ALPHA:A BETA:B\n')
                    flag = False
                if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':  # Only write over atoms
                    if line[16] != 'B' and line[26] == ' ':  # Don't allow secondary atoms
                        if line[21] == tcr_list.get('ALPHA') or line[21] == tcr_list.get('BETA'):
                            if line[21] == tcr_list.get('ALPHA') and res_alpha_count <= alpha_cut:
                                if int(line[22:26]) != previous_count_a and not flag_a:
                                    previous_count_a = int(line[22:26])
                                    flag_a = True
                                elif int(line[22:26]) != previous_count_a and flag_a:
                                    previous_count_a = int(line[22:26])
                                    res_alpha_count += 1
                                line = line[:21] + 'A' + line[22:]
                                line = line[:22] + str(res_alpha_count).rjust(4) + line[26:]
                                num = line[6:11]
                                atom_count += 1
                                if line[16] == 'A':
                                    line = line[:16] + ' ' + line[17:]
                                if res_alpha_count <= alpha_cut:
                                    output.append(line.replace(num, str(atom_count).rjust(5), 1))
                            elif line[21] == tcr_list.get('BETA') and res_beta_count <= beta_cut:
                                if int(line[22:26]) != previous_count_b and not flag_b:
                                    previous_count_b = int(line[22:26])
                                    flag_b = True
                                elif int(line[22:26]) != previous_count_b and flag_b:
                                    previous_count_b = int(line[22:26])
                                    res_beta_count += 1
                                line = line[:21] + 'B' + line[22:]
                                line = line[:22] + str(res_beta_count).rjust(4) + line[26:]
                                num = line[6:11]
                                atom_count += 1
                                if line[16] == 'A':
                                    line = line[:16] + ' ' + line[17:]
                                if res_beta_count <= beta_cut:
                                    output.append(line.replace(num, str(atom_count).rjust(5), 1))
        with open(tcr, 'w+') as f1:
            for line in output:
                f1.write(line)

    def split_mhc(self):
        """
        Creates a new PDB file with information for only the MHC of the original PDB file
        """
        tcr = self.get_pdb_id() + ".pdb"
        mhc = self.get_mhc_chain()
        helix_count = 0
        sheet_count = 0
        atom_count = 0
        compare_conect = {}
        output = []
        with open(self.file_name) as f:
            for line in f:
                left_conect = 6
                right_conect = 11
                if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                    if line[21] == mhc:
                        num = line[6:11]
                        atom_count += 1
                        output.append(line.replace(num, str(atom_count).rjust(5), 1))
                        compare_conect[num] = atom_count
                elif line[0:6] == 'HELIX ':
                    if line[19] == mhc:
                        num = line[6:10]
                        helix_count += 1
                        output.append(line.replace(num, str(helix_count).rjust(4), 1))
                elif line[0:6] == 'SHEET ':
                    if line[21] == mhc:
                        num = line[6:10]
                        sheet_count += 1
                        output.append(line.replace(num, str(sheet_count).rjust(4), 1))
                elif line[0:6] == 'HETATM':
                    if line[21] == mhc:
                        num = line[6:11]
                        atom_count += 1
                        output.append(line.replace(num, str(atom_count).rjust(5), 1))
                elif line[0:6] == 'CONECT':
                    if compare_conect.__contains__(line[left_conect:right_conect]):
                        while compare_conect.__contains__(line[left_conect:right_conect]):
                            line_update = line.replace(line[left_conect:right_conect],
                                                        str(compare_conect[line[left_conect:right_conect]]).rjust(5), 1)
                            left_conect += 5
                            right_conect += 5
                        output.append(line_update)
                else:
                    if line[0:6] != 'MASTER':
                        output.append(line)
        with open(tcr, 'w+') as f1:
            for line in output:
                f1.write(line)

    def split_p(self):
        """
        Creates a new PDB file with information for only the peptide of the original PDB file
        """
        tcr = '%s.pdb' % (self.get_pdb_id())
        peptide = self.get_peptide_chain()
        helix_count = 0
        sheet_count = 0
        atom_count = 0
        compare_conect = {}
        output = []
        with open(self.file_name) as f:
            for line in f:
                left_conect = 6
                right_conect = 11
                if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                    if line[21] == peptide:
                        num = line[6:11]
                        atom_count += 1
                        output.append(line.replace(num, str(atom_count).rjust(5), 1))
                        compare_conect[num] = atom_count
                elif line[0:6] == 'HELIX ':
                    if line[19] == peptide:
                        num = line[6:10]
                        helix_count += 1
                        output.append(line.replace(num, str(helix_count).rjust(4), 1))
                elif line[0:6] == 'SHEET ':
                    if line[21] == peptide:
                        num = line[6:10]
                        sheet_count += 1
                        output.append(line.replace(num, str(sheet_count).rjust(4), 1))
                elif line[0:6] == 'HETATM':
                    if line[21] == peptide:
                        num = line[6:11]
                        atom_count += 1
                        output.append(line.replace(num, str(atom_count).rjust(5), 1))
                elif line[0:6] == 'CONECT':
                    if compare_conect.__contains__(line[left_conect:right_conect]):
                        while compare_conect.__contains__(line[left_conect:right_conect]):
                            line_update = line.replace(line[left_conect:right_conect],
                                                       str(compare_conect[line[left_conect:right_conect]]).rjust(5), 1)
                            left_conect += 5
                            right_conect += 5
                        output.append(line_update)
                else:
                    if line[0:6] != 'MASTER':
                        output.append(line)
        with open(tcr, 'w+') as f1:
            for line in output:
                f1.write(line)

    def split_pmhc(self, update_name="..."):
        """
        Creates a new PDB file with information for only the peptide and MHC chains of the original PDB file
        Doesn't update numbering

        Parameters
        ----------
        update_name : str
            Optional naming for created PDB file
        """
        if update_name == "...":
            pmhc = 'pmhc.pdb'  # name of resulting file
        else:
            pmhc = update_name
        #peptide = self.get_peptide_chain()
        peptide = "C"
        #mhc = self.get_mhc_chain()
        mhc = "A"
        output = []
        with open(self.file_name) as f:
            for line in f:
                if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                    if line[21] == peptide or line[21] == mhc:
                        output.append(line)
                else:
                    if line[0:6] != 'MASTER':
                        output.append(line)
        with open(pmhc, 'w+') as f1:
            for line in output:
                f1.write(line)

    def split_tcr(self, update_name="...", assume_rename=False):
        """"
        Creates a new PDB file with information for only the peptide of the original PDB file
        Doesn't update numbering
        Assumes Alpha and Beta are D & E

        Parameters
        __________
        update_name : str
            Optional naming for created PDB file
        assume_rename : boolean
            Assume that the TCR chains are labeled D and E
        """
        if update_name == "...":
            tcr = 'tcr.pdb'  # name of resulting file
        else:
            tcr = update_name
        if assume_rename:  # If trimmed and renamed
            alpha = "D"
            beta = "E"
        else:
            alpha = self.get_tcr_chains()["ALPHA"]
            beta = self.get_tcr_chains()["BETA"]
        output = []
        with open(self.file_name) as f:
            for line in f:
                if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                    if line[21] == alpha or line[21] == beta:
                        output.append(line)
                else:
                    if line[0:6] != 'MASTER':
                        output.append(line)
        with open(tcr, 'w+') as f1:
            for line in output:
                f1.write(line)

    def clean_docking_count(self, rename='****'):
        """
        Provide a PDB with updated number for a TCR PDB file

        Parameters
        ----------
        rename : str
            Optional naming for created PDB file
        """
        POSSEQ = [22, 26]
        ANUM = [6, 11]
        CHAINID = 21

        if rename != '****':
            pdb_1 = rename
        else:
            pdb_1 = self.file_name  # Default naming if no input
        atom_count = 1  # Keeps track of atom number
        old_res_count = -10000
        res_count = 0  # Keeps track of residue number
        chains = []  # Chains to keep track of previous
        header = False  # Marks down header
        with open(self.file_name, "r") as i:
            with open(pdb_1, 'w+') as o:
                for line in i:
                    if line[0:6] != "MODEL " or line[0:6] != "ENDMDL":
                        if line[0:6] == 'HEADER':
                            o.write(line)
                            header = True  # Marks header
                        if header:
                            o.write('EXPDTA    DOCKING MODEL           RENUMBERED\n')
                            header = False  # Marks EXPDTA
                        if line[0:6] == 'ATOM  ':  # Only writes over atoms
                            if line[16] != 'B' and line[26] == ' ':  # Don't allow secondary atoms
                                temp_num = line[ANUM[0]:ANUM[1]]  # Saves temp old atom count
                                temp_res = int(line[POSSEQ[0]:POSSEQ[1]])  # Saves temp old res count
                                if len(chains) == 0:
                                    chains.append(line[CHAINID])
                                elif line[CHAINID] != chains[-1]:  # Catches when a new chain starts
                                    chains.append(line[CHAINID])
                                    o.write('TER\n')
                                if temp_res != old_res_count:  # Increases residue count
                                    old_res_count = int(line[POSSEQ[0]:POSSEQ[1]])
                                    res_count += 1
                                # Replaces atom count
                                line = line.replace(temp_num, str(atom_count).rjust(5), 1)
                                atom_count += 1
                                # Replaces residue count on line
                                line = line[:22] + str(res_count).rjust(4) + line[26:]
                                o.write(line)
                o.write('TER\nEND\n')

    def clean_docking_count_non_tcr(self, rename="****"):
        """
        Provide a PDB with updated number for a PDB file (Not specific to a TCR)

        Parameters
        ----------
        rename : str
            Optional naming for created PDB file
        """
        POSSEQ = [22, 26]
        ANUM = [6, 11]
        CHAINID = 21

        if rename != '****':
            renum_name = rename
        else:
            renum_name = self.file_name[:-4] + "_renum.pdb"  # Default naming if no input
        atom_count = 1  # Keeps track of atom number
        old_res_count = -10000
        res_count = 0  # Keeps track of residue number
        chains = []  # Chains to keep track of previous
        header = False  # Marks down header
        with open(self.file_name, "r") as i:
            with open(renum_name, 'w+') as o:
                for line in i:
                    if line[0:6] == 'HEADER':
                        o.write(line)
                        header = True  # Marks header
                    if header:
                        o.write('EXPDTA    DOCKING MODEL           RENUMBERED\n')
                        header = False  # Marks EXPDTA
                    if line[0:6] == 'ATOM  ':  # Only writes over atoms
                        if line[16] != 'B' and line[26] == ' ':  # Don't allow secondary atoms
                            temp_num = line[ANUM[0]:ANUM[1]]  # Saves temp old atom count
                            temp_res = int(line[POSSEQ[0]:POSSEQ[1]])  # Saves temp old res count
                            if len(chains) == 0:
                                chains.append(line[CHAINID])
                            elif line[CHAINID] != chains[-1]:  # Catches when a new chain starts
                                chains.append(line[CHAINID])
                                o.write('TER\n')
                            if temp_res != old_res_count:  # Increases residue count
                                old_res_count = int(line[POSSEQ[0]:POSSEQ[1]])
                                res_count += 1
                            # Replaces atom count
                            line = line.replace(temp_num, str(atom_count).rjust(5), 1)
                            atom_count += 1
                            # Replaces residue count on line
                            line = line[:22] + str(res_count).rjust(4) + line[26:]
                            o.write(line)
                o.write('TER\nEND\n')

    def clean_pdb(self):
        """
        Returns a reformatted PDB with just the primary TCRpMHC labeled A: MHC, B: B2M, C: Pep, D: Alpha, E: Beta
        """
        tcr_list = self.get_tcr_chains()
        atom_count = 0
        atom_flag = False  # Catch when atom and ter section is done to append tcr alpha and beta chains
        flag = False
        output = []
        other_output = []  # Hold other chains
        ab_output = []  # Temp. holds alpha beta chains until output time.
        alpha = tcr_list["ALPHA"]
        beta = tcr_list["BETA"]
        mhc = self.get_mhc_chain()
        b2m = self.get_b2m_chain()
        pep = self.get_peptide_chain()
        with open(self.file_name) as f:
            for line in f:
                if line[0:6] != 'ANISOU':  # Skip ANISOU id
                    if line[0:6] == 'HEADER':
                        output.append(line)
                        flag = True
                    if flag:
                        output.append('EXPDTA    THEORETICAL MODEL    CLEAN TCR ALPHA:D BETA:E\n')
                        flag = False
                    if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                        # Skip secondary atom positions
                        if line[16] != 'B':
                            # Selective to chains that are primary in file.
                            if line[21] == alpha or line[21] == beta\
                                    or line[21] == mhc or line[21] == b2m\
                                    or line[21] == pep:
                                if line[21] == alpha:
                                    line = line[:21] + 'D' + line[22:]
                                    ab_output.append(line)
                                elif line[21] == beta:
                                    line = line[:21] + 'E' + line[22:]
                                    ab_output.append(line)
                                elif line[21] == mhc:
                                    line = line[:21] + 'A' + line[22:]
                                    other_output.append(line)
                                elif line[21] == b2m:
                                    line = line[:21] + 'B' + line[22:]
                                    other_output.append(line)
                                elif line[21] == pep:
                                    line = line[:21] + 'C' + line[22:]
                                    other_output.append(line)
                        atom_flag = True
                    if line[0:6] != 'ATOM  ' and line[0:6] != 'TER   ' and atom_flag:
                        other_output.extend(ab_output)  # Ensures alpha and beta chains are at end.
                        for temp_line in other_output:
                            num = line[6:11]
                            atom_count += 1
                            if temp_line[16] == 'A':
                                temp_line = temp_line[:16] + ' ' + temp_line[17:]
                            output.append(temp_line.replace(num, str(atom_count).rjust(5), 1))
                        output.append("END\n")
                        break
        with open(self.file_name, 'w+') as f1:
            for line in output:
                f1.write(line)

    def fasta_TCR(self, file_name='result.fasta'):
        """
        Appends to a fasta formatted file of PDB file submitted.
        Adds sequence of alpha then beta chain

        Parameters
        ----------
        file_name : str
            Optinal naming of fasta file output
        """
        tcr_alpha_chain = self.get_amino_acid_on_chain(self.get_tcr_chains()['ALPHA'])
        tcr_beta_chain = self.get_amino_acid_on_chain(self.get_tcr_chains()['BETA'])
        total_chain = tcr_alpha_chain + tcr_beta_chain
        pdb_id = self.get_pdb_id()
        count_1 = 1
        with open(file_name, 'a+') as f:
            if tcr_alpha_chain != '' or tcr_beta_chain != '':
                f.write('>' + pdb_id + '\n')
                for aa in total_chain:
                    if count_1 % 81 != 0:
                        f.write(aa)
                        count_1 += 1
                    else:
                        f.write('\n' + aa)
                        count_1 = 2
                f.write('\n')

    def unmute_aa(self, left_aa, right_aa, chain):
        """
        Unmutes atoms of amino acid positions to untrim tcr based on left, right, and chain_id
        Assign left_aa = 0 and right_aa = 1000 for universal chain unmute

        Parameters
        ----------
        left_aa : int
            Left side of chain to mute
        right_aa : int
            Right side of chain to mute
        chain : str
        """
        range_aa = range(left_aa, right_aa + 1)
        with open(self.file_name, 'r') as r:
            data = r.readlines()
        with open(self.file_name, 'w+') as w:
            for line in data:
                if line[0:6] == 'DEATOM':
                    if line[21] == chain:
                        num_aa = int(line[22:26])
                        if num_aa in range_aa:
                            result = line.replace('DEATOM', 'ATOM  ', 1)
                            w.write(result)
                        else:
                            w.write(line)
                    else:
                        w.write(line)
                else:
                    w.write(line)

    def mute_aa(self, left_aa, right_aa, chain_id):
        """
        Mutes atoms of amino acid positions to trim tcr based on left, right, and chain_id
        Assign left_aa = 0 and right_aa = 1000 for universal chain mute

        Parameters
        ----------
        left_aa : int
            Left side of chain to mute
        right_aa : int
            Right side of chain to mute
        chain_id : str
        """
        range_aa = range(left_aa + 1, right_aa + 1)
        with open(self.file_name, 'r') as r:
            data = r.readlines()
        with open(self.file_name, 'w+') as w:
            for line in data:
                if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                    if line[21] == chain_id:
                        num_aa = int(line[22:26])
                        if num_aa in range_aa:
                            result = line.replace('ATOM  ', 'DEATOM', 1)
                            w.write(result)
                        else:
                            w.write(line)
                    else:
                        w.write(line)
                else:
                    w.write(line)

    def remove_chain(self, chain_id):
        """
         Remove chain provided based on ID
        Parameters
        ----------
        chain_id : str
        """
        with open(self.file_name, 'r') as r:
            data = r.readlines()
        with open(self.file_name, 'w+') as w:
            for line in data:
                if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                    if line[21] != chain_id.upper():
                        w.write(line)

    def trim_chain(self, chain_id, cutoff):
        """
        Trims chain submitted with cutoff provided of AA count

        Parameters
        ----------
        chain_id : str
        cutoff : int
            Position in chain to cut
        """
        output = []
        with open(self.file_name, "r") as i:
            for line in i:
                if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':  # Only write over atoms
                    if line[16] != 'B' and line[26] == ' ':  # Don't allow secondary atoms
                        if line[21] == chain_id.upper() and int(line[22:26]) <= cutoff or line[21] != chain_id.upper():
                            output.append(line)
                else:
                    output.append(line)
        with open(self.file_name, "w") as o:
            for line in output:
                o.write(line)

    def split_chains(self, chains_in, suffix, dir_location='****'):
        """
        Splits chains given and renames with given suffix

        Parameters
        ----------
        chains_in : str
        suffix : str
            Suffix of created file
        dir_location : str
            Optional submission of location of PDB file
        """
        if dir_location == '****':
            new_pdb = self.get_file_name()
        else:
            new_pdb = dir_location + self.get_pdb_id() + suffix + ".pdb"
        output = []
        with open(self.file_name) as f:
            for line in f:
                if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                    if line[21] in chains_in.upper():
                        output.append(line)
                else:
                    if line[0:6] != 'MASTER':
                        output.append(line)
        with open(new_pdb, 'w+') as f1:
            for line in output:
                f1.write(line)

    def superimpose(self, ref_pdb, target_order, ref_order, new_name_in="..."):
        """
        Superimpose two PDBs
        Read in target and reference structure, superimpose target to reference, save superimposed target structure

        Parameters
        ----------
        ref_pdb : str
            Location of reference PDB
        target_order : str
            Order of chains to compare in superimposing structures for target
        ref_order : str
            Order of chains to compare in superimposing structures for reference
        new_name_in : str
            Optional name in for resulting superimposed positioning of target structure

        Returns
        -------
            super_imposer.rms: float
                RMSD value for all-atom RMSD
        """
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.gap_score = -100.00
        aligner.match_score = 2.0
        aligner.mismatch_score = 0.0

        # Collect target pdb seq
        target_pdb = self.get_file_name()
        target_seq = {}  # Dictionary chain_id:seq
        for chain in self.get_chains():
            target_seq[chain] = self.get_amino_acid_on_chain(chain)

        # Collect reference pdb seq
        self.set_file_name(ref_pdb)
        ref_seq = {}  # Dictionary chain_id:seq
        for chain in self.get_chains():
            ref_seq[chain] = self.get_amino_acid_on_chain(chain)
        # Set back to target pdb
        self.set_file_name(target_pdb)

        # Starting positions of each AA from target and reference
        # 'target' -> 'chain' -> [start, end]
        start_pos = {"target": {}, "reference": {}}
        for chain in list(target_order):
            # print(chain)
            start_pos['target'][chain] = []
        for chain in list(ref_order):
            # print(chain)
            start_pos['reference'][chain] = []

        # Produce alignments
        target_pos = list(target_order)
        ref_pos = list(ref_order)
        for chain_pos in range(0, len(target_pos)):
            # Align with chain in target_order and ref_order
            target_chain = target_seq[target_pos[chain_pos]]
            reference_chain = ref_seq[ref_pos[chain_pos]]
            temp_align = aligner.align(target_chain, reference_chain)
            # Uncomment these lines to show alignments
            # print("Target Chain: " + target_pos[chain_pos])
            # print("Reference Chain: " + ref_pos[chain_pos])
            # print(temp_align[0])
            start_pos['target'][target_pos[chain_pos]] = [temp_align[0].path[0][0], temp_align[0].path[1][0]]
            start_pos['reference'][ref_pos[chain_pos]] = [temp_align[0].path[0][1], temp_align[0].path[1][1]]

        # Initialize parser
        parser = Bio.PDB.PDBParser(QUIET=True)

        # Gather Structures
        ref_structure = parser.get_structure("reference", ref_pdb)
        target_structure = parser.get_structure("target", self.file_name)

        # Collect structures
        ref_model = ref_structure[0]
        target_model = target_structure[0]

        # Collecting alpha carbons
        ref_atoms = []
        target_atoms = []
        # Reference Model
        for ref_chain in ref_model:
            chain = ref_chain.__repr__().split("=")[1].split(">")[0]
            # print(ref_chain)
            # print(chain)
            if chain in ref_pos:
                first_in_chain = True
                skip_switch = False
                start = 0
                end = 0
                last_count = 0
                for ref_res in ref_chain:
                    if ref_res.get_resname() != "HOH" and 'CA' in ref_res:
                        if first_in_chain:
                            offset = ref_res.get_id()[1]
                            start = start_pos['reference'][chain][0] + offset
                            end = start_pos['reference'][chain][1] + offset
                            first_in_chain = False
                        if ref_res.get_id()[1] in range(start, end):
                            if skip_switch:
                                if ref_res.get_id()[1] != last_count + 1:  # Catch when count skips a number
                                    end += ref_res.get_id()[1] - last_count - 1
                            ref_atoms.append(ref_res['CA'])
                            last_count = ref_res.get_id()[1]
                            skip_switch = True
        # Target Model
        for target_chain in target_model:
            chain = target_chain.__repr__().split("=")[1].split(">")[0]
            if chain in target_pos:
                first_in_chain = True
                skip_switch = False
                start = 0
                end = 0
                last_count = 0
                for target_res in target_chain:
                    if target_res.get_resname() != "HOH" and 'CA' in target_res:
                        if first_in_chain:  # If chain doesn't start count at position 0
                            offset = target_res.get_id()[1]
                            start = start_pos['target'][chain][0] + offset
                            end = start_pos['target'][chain][1] + offset
                            first_in_chain = False
                        if target_res.get_id()[1] in range(start, end):
                            if skip_switch:
                                if target_res.get_id()[1] != last_count + 1:  # Catch when count skips a number
                                    end += target_res.get_id()[1] - last_count - 1
                            target_atoms.append(target_res['CA'])
                            last_count = target_res.get_id()[1]
                            skip_switch = True
        # print(len(ref_atoms))
        # print(len(target_atoms))
        # Truncate long list
        if len(ref_atoms) > len(target_atoms):
            ref_atoms = ref_atoms[:len(target_atoms)]
        elif len(target_atoms) > len(ref_atoms):
            target_atoms = target_atoms[:len(ref_atoms)]
        # Superimposing
        super_imposer = Bio.PDB.Superimposer()
        super_imposer.set_atoms(ref_atoms, target_atoms)
        super_imposer.apply(target_model.get_atoms())

        # Save structure
        io = Bio.PDB.PDBIO()
        if __name__ == "__main__":  # Dont run if being used as package
            print("RMSD: " + str(super_imposer.rms) + " Atoms Pulled: " + str(len(target_atoms)))
        io.set_structure(target_structure)
        if new_name_in != "...":
            new_name = new_name_in
        else:
            new_name = self.get_file_name().split(".")[0] + "_aligned.pdb"
        io.save(new_name)
        return super_imposer.rms

    def rmsd(self, ref_pdb, target_order, ref_order, ca=False, mute=False):
        """
        Calculate RMSD values between two PDBs - based on aligned aa's. Optional Carbon Alpha RMSD as well.

        Parameters
        ----------
        ref_pdb : str
            Reference PDB - can be same PDB as target
        target_order : str
            Target chains to consider for RMSD calculation - must be in same order
        ref_order : str
            Reference chains to consider for RMSD calculation - must be in same order
        ca : boolean
            True - only Alpha Carbon RMSD ; False - all-atom RMSD
        mute : boolean
            Don't print text, only return when called

        Returns
        -------
        rmsd : float
        """
        target_atoms = {}  # chain: atoms
        ref_atoms = {}  # chain: atoms
        target_aa = {}  # chain: aa's
        ref_aa = {}  # chain: aa's
        rmsd_array = {"target": [], "ref": []}  # List of cords in order
        current_pdb = self.get_file_name()
        # Collect target atoms
        for chain in target_order:
            # Collect atoms - contains xyz and aa
            target_atoms[chain] = self.get_atoms_on_chain(chain)
            # Collect just residue list
            target_aa[chain] = self.get_amino_acid_on_chain(chain)
        # Collect reference atoms
        self.set_file_name(ref_pdb)
        for chain in ref_order:
            # Collect atoms - contains xyz and aa
            ref_atoms[chain] = self.get_atoms_on_chain(chain)
            # Collect just residue list
            ref_aa[chain] = self.get_amino_acid_on_chain(chain)
        # Return to previous PDB
        self.set_file_name(current_pdb)
        # Run alignment
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.gap_score = -100.00
        aligner.match_score = 2.0
        aligner.mismatch_score = 0.0
        # Loop through each paired chain
        for position in range(len(target_order)):
            # run alignment
            temp_align = aligner.align(target_aa[target_order[position]], ref_aa[ref_order[position]])
            # Collect info needed - start and end positions of alignments ex. [0, 101]
            target_chain_info = [temp_align[0].path[0][0], temp_align[0].path[1][0]]
            ref_chain_info = [temp_align[0].path[0][1], temp_align[0].path[1][1]]
            # Collect XYZ cords. for target
            last_pos = 0
            count = -1
            for atom in target_atoms[target_order[position]]:
                if atom["atom_id"] == "CA" or not ca:  # Check if not carbon alpha
                    if atom["atom_id"][0] != "H":  # Skip hydrogen atoms added while docking
                        if count == -1:
                            last_pos = atom["comp_num"]
                            count = 0
                        if atom["comp_num"] != last_pos:
                            count += 1
                        if count in range(target_chain_info[0], target_chain_info[1]):
                            rmsd_array["target"].append([atom["X"], atom["Y"], atom["Z"]])
            # Collect XYZ cords. for ref
            count = -1
            for atom in ref_atoms[ref_order[position]]:
                if atom["atom_id"] == "CA" or not ca:  # Check if not Carbon alpha
                    if atom["atom_id"][0] != "H":  # Skip hydrogen atoms added while docking
                        if count == -1:
                            last_pos = atom["comp_num"]
                            count = 0
                        if atom["comp_num"] != last_pos:
                            count += 1
                        if count in range(ref_chain_info[0], ref_chain_info[1]):
                            rmsd_array["ref"].append([atom["X"], atom["Y"], atom["Z"]])
        # print(rmsd_array)
        # print(len(rmsd_array["target"]))
        # print(len(rmsd_array["ref"]))
        # Calculate RMSD
        sum_ = 0
        for position in range(len(rmsd_array["target"])):
            tar = rmsd_array["target"][position]
            ref = rmsd_array["ref"][position]
            sum_ += (pow(tar[0] - ref[0], 2) + pow(tar[1] - ref[1], 2) + pow(tar[2] - ref[2], 2))
        rmsd = sqrt(sum_/len(rmsd_array["target"]))
        if not mute:
            print("RMSD: " + str(rmsd))
        return rmsd

    def get_center(self):
        # Pull atoms from each chain: default_atoms {'CHAIN_ID': [list of atoms dictionaries]}
        default_atoms = {}
        cords_dic = {'X': [], 'Y': [], 'Z': []}
        for chain in self.get_chains():
            default_atoms[chain] = self.get_atoms_on_chain(chain)
            for atom in default_atoms[chain]:
                for key in cords_dic.keys():
                    cords_dic[key].append(atom[key])
        # Calculate center
        centroid = [round(statistics.mean(cords_dic['X']), 3), round(statistics.mean(cords_dic['Y']), 3),
                    round(statistics.mean(cords_dic['Z']), 3)]
        return centroid

    def center(self, new_name_in="..."):
        """
        Gather all atoms in the file and calculate the average XYZ coord to find the center of the structure

        Parameters
        ----------
        new_name_in : str
            Optional name for new PDB created by method with centered structure
        """
        atoms = []
        full_atom = []
        # Collect atom information from PDB
        for chain in self.get_chains():
            chain_atoms = self.get_atoms_on_chain(chain)
            for atom in chain_atoms:
                # Append atom_line with full atom information
                full_atom.append(atom)
                # Append atoms with just the coordinates of each atom
                atoms.append([atom['X'], atom['Y'], atom['Z']])
        # Convert to array
        atom_array = np.array(atoms)
        # Calculate avg for each axis
        current_avg = np.mean(atom_array, axis=0)
        # Center the coordinates
        new_array = []
        for set_pos in atom_array:
            x = set_pos - current_avg
            new_array.append(x)
        # Calculate PCA
        pca = PCA(n_components=3)
        pca.fit(new_array)
        transpose = np.transpose(pca.components_)
        # Determinate of components matrix - Corrects if determinate is negative
        determinate = np.linalg.det(transpose)
        if determinate < 0:
            for position in transpose:
                position[0] = position[0] * -1
        # Rotate x-axis
        # Always have N-terminus in positive coordinates (Fixes flips on x-axis)
        test_x = np.matmul(new_array[-1], transpose)
        if test_x[1] < 0:
            r = Rotation.from_euler('x', 180, degrees=True)
            transpose = r.apply(transpose)
        # Rotate y-axis
        # Always have first chain on left side | alpha on left side (Fixes flips on y-axis)
        # Determines by looking at last atom which should be on the Beta chain
        test_y = np.matmul(new_array[-1], transpose)
        if test_y[0] < 0:
            r = Rotation.from_euler('y', 180, degrees=True)
            transpose = r.apply(transpose)
        # Multiply by EV
        new_cords = np.matmul(new_array, transpose)
        # Replace XYZ coordinates
        axis = ['X', 'Y', 'Z']
        for num in range(0, len(full_atom)):
            for position in range(0, len(axis)):
                full_atom[num][axis[position]] = new_cords[num][position]
        # Write to new file
        if new_name_in != "...":
            new_name = new_name_in
        else:
            new_name = self.get_file_name().split("/")[-1].split(".")[0] + "_center.pdb"
        with open(new_name, "w") as f:
            # Send to reconstruct atom lines
            f.write(self.rebuild_atom_line(full_atom))

    def join(self, pdb_1, pdb_2, new_name):
        """
        Joins together two PDB files by appending first PDBs atoms to second PDBs atoms

        Parameters
        ----------
        pdb_1 : str
            Location and name of PDB 1
        pdb_2 : str
            Location and name of PDB 2
        new_name : str
            Name of new file
        """
        atoms_lines = []
        pdbs = [pdb_1, pdb_2]
        for pdb in pdbs:
            with open(pdb, "r") as f1:
                for line in f1:
                    if line[0:6] == "ATOM  " or line[0:6] == "TER   ":
                        atoms_lines.append(line)
        with open(new_name, "w") as f2:
            for line in atoms_lines:
                f2.write(line)
        return new_name

    def reorder_chains(self, chain_order):
        """
        Update the chain order. Must send in a list with identical number of chains

        Parameters
        ----------
        chain_order : str
            Order in which you want the chains to be in the file
        """
        current = self.get_chains()
        chain_info = {}
        new_order = []
        for chain in current:
            chain_info[chain] = self.get_atoms_on_chain(chain)
        for chain in list(chain_order):
            for atom in chain_info[chain]:
                new_order.append(atom)
        with open(self.file_name, "w") as f1:
            f1.write(self.rebuild_atom_line(new_order))

    def update_label(self, label_dic):
        """
        Update the labels based on submitted chain dictionary
        Ex of dictionary: {'A':'D', 'B':'E'}  A gets replaced with D and B gets replaced with E

        Parameters
        ----------
        label_dic : dict
            Dictionary of label names that need to be adjusted
        """
        current = self.get_chains()
        chain_info = {}
        new_order = []
        for chain in current:
            chain_info[chain] = self.get_atoms_on_chain(chain)
        for chain in chain_info:
            for atom in chain_info[chain]:
                atom['chain_id'] = label_dic[chain]
                new_order.append(atom)
        with open(self.file_name, 'w') as f1:
            f1.write(self.rebuild_atom_line(new_order))

    # Below CDR methods are adapted from Ryan Ehrlich's code
    def pull_cdr(self):
        """
        Primary method
        Returns CDR1a, CDR2a, CDR2.5a, CDR3, CDR1b, CDR2b, CDR2.5b, CDR3b

        Returns
        -------
        alpInds
        betInds
        """
        from getCDRs import cdr_loops
        trav, trbv = cdr_loops()  # get germline dictionary
        # loops = CdrLoopInfo(self.get_file_name())
        tcrDict = self.getLines()  # return residues identity and positions
        aseq, bseq, apos, bpos = self.getSeqs(tcrDict)  # return pdb seqs and positions
        # CDR loop identity and residue numbering
        alpInds, betInds = self.loopPositions(aseq, apos, trav), self.loopPositions(bseq, bpos, trbv)
        return alpInds, betInds

    def getLines(self):
        """
        Helper method to pull_cdr()

        Returns
        -------
        seq_dict
        """
        alpha = ['D']
        beta = ['E']
        data = open(self.get_file_name(), 'r')
        seq_dict = {'alpha': [], 'beta': []}
        for line in data:
            if line[:4] == "ATOM":
                if line[21] in alpha:
                    resID = line[POSRES[0]:POSRES[1]] + "_" + line[POSSEQ[0]:POSSEQ[1]] + "_" + line[POSCHAIN]
                    resID = resID.replace(" ", "")
                    if resID not in seq_dict['alpha']:
                        seq_dict['alpha'].append(resID)
                if line[21] in beta:
                    resID = line[POSRES[0]:POSRES[1]] + "_" + line[POSSEQ[0]:POSSEQ[1]] + "_" + line[POSCHAIN]
                    resID = resID.replace(" ", "")
                    if resID not in seq_dict['beta']:
                        seq_dict['beta'].append(resID)
        return seq_dict

    def getSeqs(self, dict):
        """
        Helper method to pull_cdr()
        return seq of pdb chains (tcr) as well as pos_chain

        Parameters
        ----------
        dict : dict
        """
        # triple letter abbreviation for amino acids
        aa_names = [
            "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
            "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
        ]

        # single letter abbreviation for amino acids
        aa_singleletter = "ACDEFGHIKLMNPQRSTVWY"
        aa_tripleletter = aa_names
        keys = dict.keys()
        alp, bet, alp_inds, bet_inds = [], [], [], []
        for key in keys:
            seq = dict.get(key)
            for r in range(len(seq)):
                res = seq[r][:3]
                ind = aa_tripleletter.index(res)
                res_ = aa_singleletter[ind]
                if key == 'alpha':
                    alp.append(res_)
                    alp_inds.append(seq[r][4:-2])
                if key == 'beta':
                    bet.append(res_)
                    bet_inds.append(seq[r][4:-2])
        return ''.join(alp), ''.join(bet), alp_inds, bet_inds

    def getCDR3(self, seq, germline):
        """
        Helper method to pull_cdr()

        Parameters
        ----------
        seq : str
        germline : str

        Returns
        -------
        cdr3 : str
        """
        start = germline.rfind('C')  # get ind of C from CDR dict
        check = germline.find(seq[5:10])  # reference start point in germline
        diff = check - 5  # difference between seq and germline start
        seq_ = seq[start - diff:]  # adjusted seq to match gerline
        rev = []
        for s in reversed(seq_):  # make reverse seq_
            rev.append(s)
        rev = ''.join(rev)
        end = 0
        for i in range(len(rev) - 2):  # find end of CDR3 using G*G pattern starting from the end of seq_ (rev)
            if rev[i] == 'G':
                if rev[i + 2] == 'G':
                    end = i
                    break
        if end != 0:  # if end does not equal 0, trim using G*G pattern
            end2 = len(seq_) - end  # end of G*G patter with CDR3
            cdr3 = seq_[:end2 - 3]  # trim G*G pattern
        if rev[0] == 'F' or end == 0:  # if end = 0, end of seq_ is the end of CDR3
            cdr3 = seq_
        return (cdr3)

    def loopPositions(self, seq, resNums, trv):
        """
        Helper method to pull_cdr()

        Parameters
        ----------
        seq
        resNums
        trv

        Returns
        -------
        loopInfo
        """
        cdrList = []
        for key in trv:
            loops = trv.get(key)
            check = [''.join(loops[0]), ''.join(loops[1]), ''.join(loops[2])]
            if all(x in seq for x in check):  # if all germline CDRs in seq, continue
                cdr3 = self.getCDR3(seq, ''.join(loops[4]))  # get CDR3
                cdrList = [''.join(loops[0]), ''.join(loops[1]), ''.join(loops[2]), cdr3]
                break
        loopInfo = []  # cdr loop with res number range (corresponds to pdb)
        for loop in cdrList:  # get loops and loop res number ranges
            ind = seq.find(loop)
            start, end, = resNums[ind], resNums[ind + len(loop) - 1]
            loopInfo.append([loop, int(start), int(end)])
        return (loopInfo)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="Full crystal structure", type=str)
    parser.add_argument("--get_tcr", help="Returns the TCR chain ids", default=False, action="store_true")
    parser.add_argument("--get_mhc", help="Returns the mhc chain id", default=False, action="store_true")
    parser.add_argument("--renum", help="Updates numbering of TCR for docking", default=False, action="store_true")
    parser.add_argument("--renum2", help="Updates numbering for docking with non-tcr files", default=False,
                        action="store_true")
    parser.add_argument("--trim", help="Trim TCR chains to only contain variable region", default=False,
                        action="store_true")
    parser.add_argument("--mhc_split", help="Only provide MHC chains", default=False, action="store_true")
    parser.add_argument("--peptide_split", help="Only provide MHC chains", default=False, action="store_true")
    parser.add_argument("--clean_tcr_split", help="Only provide TCR chains", default=False, action="store_true")
    parser.add_argument("--tcr_split", help="Only provides TCR chains", default=False, action="store_true")
    parser.add_argument("--tcr_split_default", help="Only provides TCR chains, assumes DE", default=False,
                        action="store_true")
    parser.add_argument("--pmhc_split", help="Only provides pMHC chains", default=False, action="store_true")
    parser.add_argument("--peptide", help="Get peptide chain", default=False, action="store_true")
    parser.add_argument("--mhc", help="Get mhc chain", default=False, action="store_true")
    parser.add_argument("--pmhc", help="Get pmhc chains", default=False, action="store_true")
    parser.add_argument("--alpha", help="Get alpha chain", default=False, action="store_true")
    parser.add_argument("--beta", help="Get beta chain", default=False, action="store_true")
    parser.add_argument("--resolution", help="Get resolution", default=False, action="store_true")
    parser.add_argument("--clean_pdb", help="Updated to updated labeling and chain order", default=False,
                        action="store_true")
    parser.add_argument("--align", help="(Align) Superimpose this reference structure to submitted pdb", type=str)
    parser.add_argument("--rmsd", help="(rmsd) Calculate RMSD of all-atoms or Carbon alpha", type=str)
    parser.add_argument("--carbon", help="(rmsd) Change to carbon alpha RMSD calculation", action="store_true",
                        default=False)
    parser.add_argument("--tar_chains", help="(Align|rmsd) Chains from target to match with reference", type=str)
    parser.add_argument("--ref_chains", help="(Align|rmsd) Chains from reference to match with target", type=str)
    parser.add_argument("--center", help="Center TCR to cord. 0,0,0", action="store_true", default=False)
    parser.add_argument("--reorder", help="Reorder chains based on string provided (case sensitive)", type=str)
    parser.add_argument("--pull_cdr", help="Pull CDRs", action="store_true", default=False)
    return parser.parse_args()


####################
#     Controls     #
####################
def main():
    args = parse_args()
    pdb = PdbTools3(args.pdb)
    if args.get_tcr:
        print(pdb.get_tcr_chains())
    if args.get_mhc:
        print(pdb.get_mhc_chain())
    if args.mhc_split:
        pdb.split_mhc()
    if args.trim:
        pdb.clean_tcr_count_trim()
    if args.clean_tcr_split:
        pdb.clean_tcr()
    if args.peptide_split:
        pdb.split_p()
    if args.tcr_split:
        pdb.split_tcr()
    if args.tcr_split_default:
        pdb.split_tcr("...", True)
    if args.pmhc_split:
        pdb.split_pmhc()
    if args.renum:
        pdb.clean_docking_count()
    if args.renum2:
        pdb.clean_docking_count_non_tcr()
    if args.peptide:
        print(pdb.get_peptide_chain())
    if args.mhc:
        print(pdb.get_mhc_chain())
    if args.alpha:
        print(pdb.get_tcr_chains()['ALPHA'])
    if args.beta:
        print(pdb.get_tcr_chains()['BETA'])
    if args.resolution:
        print(pdb.get_resolution())
    if args.clean_pdb:
        pdb.clean_pdb()
    if args.align:
        pdb.superimpose(args.align, args.tar_chains, args.ref_chains)
    if args.rmsd:
        pdb.rmsd(args.rmsd, args.tar_chains, args.ref_chains, args.carbon)
    if args.center:
        if os.path.isdir(args.pdb):
            os.mkdir("Results")
            for each in os.listdir(args.pdb):
                if each.endswith(".pdb"):
                    print(each.split(".")[0])
                    pdb.set_file_name(args.pdb + "/" + each)
                    pdb.center("Results/" + each.split(".")[0] + "_center.pdb")
        else:
            pdb.center()
    if args.reorder:
        pdb.reorder_chains(args.reorder)
    if args.pull_cdr:
        # from getCDRs import cdr_loops
        print(pdb.pull_cdr())


if __name__ == '__main__':
    main()
