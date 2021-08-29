# This file is a part of the TRain program
# Author: Austin Seamann & Dario Ghersi
# Version: 1.00
# Last Updated: August 4th, 2021
from math import sqrt
from Bio import pairwise2
import warnings
from Bio import BiopythonDeprecationWarning
with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
    from Bio.SubsMat import MatrixInfo as matlist
import argparse


class PdbTools3:
    # initialize PdbTools
    def __init__(self, file="..."):
        self.file_name = file
        self.test_list = {}

    # method for changing PDB file
    def set_file_name(self, file_name_in):
        self.file_name = file_name_in

    # Returns file name currently in use
    def get_file_name(self):
        return self.file_name

    # Returns the pdbID of file
    def get_pdb_id(self):
        with open(self.file_name, 'r') as file:
            output = file.readline()
            return output[62:66].lower()

    # Returns a list of all chains in PDB file
    def get_chains(self):
        chains = []
        with open(self.file_name, 'r') as file:
            for line in file:
                if line[0:6] == 'ATOM  ':
                    if not chains.__contains__(line[21]):
                        chains.append(line[21])
        return chains

    def get_resolution(self):
        value = ''
        output = -1
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

    # Returns a string of amino acids in a specific chain as a string in single letter notation
    def get_amino_acid_on_chain(self, chain):
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

    # Returns a dictionary for the first atom of a chain.
    def first_atom_on_chain(self, chain):
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

    # Returns a dictionary with elements related to specific atom number entered
    def get_atom(self, atom_num):
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

    # Returns the Euclidean distance between two atoms based with atom_id being sent in as the parameter
    def euclidean_of_atoms(self, atom_num_1, atom_num_2):
        atom_1 = self.get_atom(atom_num_1)
        atom_2 = self.get_atom(atom_num_2)
        euclidean_distance = sqrt((atom_2['X'] - atom_1['X'])**2 + (atom_2['Y'] - atom_1['Y'])**2
                                  + (atom_2['Z'] - atom_1['Z'])**2)
        return euclidean_distance

    def get_atoms_on_chain(self, chain):
        atoms = []
        with open(self.file_name, 'r') as file:
            for line in file:
                if line[0:6] == 'ATOM  ':
                    if line[21] == chain.upper() and len(line) >= 76:
                        atoms.append({'atom_num': int(line[6:11]), 'atom_id': line[13:16].strip(),
                                'atom_comp_id': line[17:20],
                                'chain_id': line[21], 'comp_num': int(line[22:26]), 'X': float(line[31:38]),
                                'Y': float(line[38:46]), 'Z': float(line[46:54]), 'occupancy': float(line[55:60]),
                                'B_iso_or_equiv': float(line[60:66]), 'atom_type': line[77]})
                    elif line[21] == chain.upper() and len(line) >= 76:
                        atoms.append({'atom_num': int(line[6:11]), 'atom_id': line[13:16].strip(),
                                'atom_comp_id': line[17:20],
                                'chain_id': line[21], 'comp_num': int(line[22:26]), 'X': float(line[31:38]),
                                'Y': float(line[38:46]), 'Z': float(line[46:54]), 'occupancy': float(line[55:60])})
        return atoms

    # Returns alpha and beta chain IDs based on seq. alignment to 1a07 PDB entry chains. Confirms that it is a
    # partnering TCR chain
    def get_tcr_chains(self):
        matrix = matlist.blosum62
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
        for chain in chains:
            score_alpha = pairwise2.align.globaldx(self.get_amino_acid_on_chain(chain), alpha_chain[0], matrix,
                                                   score_only=True, penalize_end_gaps=(False, False))
            tmp_alpha.append([float(score_alpha), chain])
            score_beta = pairwise2.align.globaldx(self.get_amino_acid_on_chain(chain), beta_chain[0], matrix,
                                                  score_only=True, penalize_end_gaps=(False, False))
            tmp_beta.append([float(score_beta), chain])
        alpha = sorted(tmp_alpha)
        beta = sorted(tmp_beta)
        result['ALPHA'] = alpha[-1][1]
        result['BETA'] = beta[-1][1]
        atom = self.first_atom_on_chain(result['ALPHA'])
        position = -1
        while True:
            atom2 = self.first_atom_on_chain(result['BETA'])
            # Distance between 1st atom in each chain
            distance1 = self.euclidean_of_atoms(atom['atom_num'], atom2['atom_num'])
            # Distance between 1st and 125 atoms in each chain
            distance2 = self.euclidean_of_atoms(atom['atom_num'], atom2['atom_num'] + 125)
            # Determines if TCR chains are within typical distance.
            if distance1 <= 46 and abs(distance1 - distance2) <= 15:
                self.test_list[self.get_file_name()] = abs(distance1 - distance2)
                result['BETA'] = beta[position][1]
                break
            elif position == (len(beta) * -1):
                print(self.get_file_name() + ' : Possibly does not contain both TCR chains')
                result['BETA'] = beta[-1][1]
                break
            else:
                position -= 1
                result['BETA'] = beta[position][1]
        return result

    # Returns the amino acid sequence of either 'ALPHA' or 'BETA' chain as single letter AA abbreviation
    def get_tcr_amino_seq(self, tcr_type_in):
        tcr_dict = self.get_tcr_chains()
        for key in tcr_dict:
            if key == tcr_type_in:
                return self.get_amino_acid_on_chain(tcr_dict[key])

    # Returns the atoms on the peptide, peptide is determined by smallest chain
    def get_atoms_on_peptide(self):
        atoms = []
        self.set_record_type('ATOM')
        file_atoms = self.record_report().splitlines()
        chain = self.get_peptide_chain()
        for line in file_atoms:
            if line[15] == chain:
                atoms.append(self.get_atom(int(line[0:5])))
        atoms.pop(len(atoms) - 1)
        return atoms

    # Returns the AA chain that is the peptide of the pMHC complex, based on the smallest chain in file
    def get_peptide_chain(self):
        chains = self.get_chains()
        chain_dic = {}
        length = 10000
        peptide = ''
        for chain in chains:
            chain_dic[chain] = self.get_amino_acid_on_chain(chain)
        for chain in chain_dic:
            if len(chain_dic[chain]) < length:
                length = len(chain_dic[chain])
                peptide = chain
        return peptide

    # Returns the MHC chain - based on COMPND section labeling HISTOCOMPATIBILITY ANTIGEN
    # WILL LATER CONVERT THIS TO HOW I CHECK FOR TCR CHAINS
    def get_mhc_chain(self):
        matrix = matlist.blosum62
        # Hard coded mhc chain
        mhc_chain = 'GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTV'\
                    'QRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMT'\
                    'HHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWE'
        chains = self.get_chains()
        tmp_mhc = []
        for chain in chains:
            score_mhc = pairwise2.align.globaldx(
                self.get_amino_acid_on_chain(chain), mhc_chain, matrix,
                score_only=True, penalize_end_gaps=(False, False))
            tmp_mhc.append([float(score_mhc), chain])
        mhc = sorted(tmp_mhc)
        return mhc[-1][1]

    def renumber_docking(self, rename="****"):
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

    # Returns a reformatted PDB with just the TCR atom cord. and has relabeled chains with ALPHA = A; BETA = B
    def clean_tcr(self, dir_start='****'):
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

    # Returns a reformatted PDB with just the TCR atom cord. and has relabeled chains with ALPHA = A; BETA = B
    # Also uses trimmed TCR seq.
    # New count on amino acids so each chain starts at 1
    def clean_tcr_count_trim(self, dir_start='****'):
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

    # Creates a new PDB file with information for only the MHC of the original PDB file
    def split_mhc(self):
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

    # Creates a new PDB file with information for only the peptide of the original PDB file
    def split_p(self):
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

    # Creates a new PDB file with information for only the peptide of the original PDB file
    # Doesn't update numbering
    def split_pmhc(self, update_name="..."):
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

    # Creates a new PDB file with information for only the peptide of the original PDB file
    # Doesn't update numbering
    # Assumes Alpha and Beta are D & E
    def split_tcr(self, update_name="..."):
        if update_name == "...":
            tcr = 'tcr.pdb'  # name of resulting file
        else:
            tcr = update_name
        #alpha = self.get_tcr_chains()["ALPHA"]
        #beta = self.get_tcr_chains()["BETA"]
        alpha = "D"
        beta = "E"
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
        POSSEQ = [22, 26]
        POSCHAIN = 21
        ANUM = [6, 11]
        CHAINID = 21

        if rename != '****':
            pdb_1 = rename
        else:
            pdb_1 = self.pdb  # Default naming if no input
        atom_count = 1  # Keeps track of atom number
        old_res_count = -10000
        res_count = 0  # Keeps track of residue number
        chains = []  # Chains to keep track of previous
        header = False  # Marks down header
        with open(self.pdb, "r") as i:
            with open(pdb_1, 'w+') as o:
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

    # Returns a reformatted PDB with just the primary TCRpMHC labeled A: MHC, B: B2M, C: Pep, D: Alpha, E: Beta
    def clean_pdb(self):
        tcr_list = self.get_tcr_chains()
        atom_count = 0
        flag = False
        output = []
        alpha = tcr_list["ALPHA"]
        beta = tcr_list["BETA"]
        mhc = self.get_mhc_chain()
        b2m = self.get_b2m_chain()
        pep = self.get_peptide_chain()
        with open(self.file_name) as f:
            for line in f:
                if line[0:6] == 'HEADER':
                    output.append(line)
                    flag = True
                if flag:
                    output.append('EXPDTA    THEORETICAL MODEL    CLEAN TCR ALPHA:D BETA:E\n')
                    flag = False
                if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                    if line[16] != 'B':
                        # Selective to chains that are primary in file.
                        if line[21] == alpha or line[21] == beta\
                                or line[21] == mhc or line[21] == b2m\
                                or line[21] == pep:
                            if line[21] == alpha:
                                line = line[:21] + 'D' + line[22:]
                            elif line[21] == beta:
                                line = line[:21] + 'E' + line[22:]
                            elif line[21] == mhc:
                                line = line[:21] + 'A' + line[22:]
                            elif line[21] == b2m:
                                line = line[:21] + 'B' + line[22:]
                            elif line[21] == pep:
                                line = line[:21] + 'C' + line[22:]
                            num = line[6:11]
                            atom_count += 1
                            if line[16] == 'A':
                                line = line[:16] + ' ' + line[17:]
                            output.append(line.replace(num, str(atom_count).rjust(5), 1))
        with open(self.file_name, 'w+') as f1:
            for line in output:
                f1.write(line)

    # Appends to a fasta formatted file of PDB file submitted.
    # Adds sequence of alpha then beta chain
    def fasta_TCR(self, file_name='result.fasta'):
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

    # Unmutes atoms of amino acid positions to untrim tcr based on left, right, and chain_id
    # Assign left_aa = 0 and right_aa = 1000 for universal chain unmute
    def unmute_aa(self, left_aa, right_aa, chain):
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

    # Assign left_aa = 0 and right_aa = 1000 for universal chain mute
    def mute_aa(self, left_aa, right_aa, chain_id):
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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="Full crystal structure", type=str)
    parser.add_argument("--renum", help="Updates numbering of TCR for docking", default=False, action="store_true")
    parser.add_argument("--trim", help="Trim TCR chains to only contain variable region", default=False,
                        action="store_true")
    parser.add_argument("--mhc_split", help="Only provide MHC chains", default=False, action="store_true")
    parser.add_argument("--peptide_split", help="Only provide MHC chains", default=False, action="store_true")
    parser.add_argument("--clean_tcr_split", help="Only provide TCR chains", default=False, action="store_true")
    parser.add_argument("--tcr_split", help="Only provides TCR chains", default=False, action="store_true")
    parser.add_argument("--pmhc_split", help="Only provides pMHC chains", default=False, action="store_true")
    parser.add_argument("--peptide", help="Get peptide chain", default=False, action="store_true")
    parser.add_argument("--mhc", help="Get mhc chain", default=False, action="store_true")
    parser.add_argument("--pmhc", help="Get pmhc chains", default=False, action="store_true")
    parser.add_argument("--alpha", help="Get alpha chain", default=False, action="store_true")
    parser.add_argument("--beta", help="Get beta chain", default=False, action="store_true")
    parser.add_argument("--resolution", help="Get resolution", default=False, action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    pdb = PdbTools3(args.pdb)
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
    if args.pmhc_split:
        pdb.split_pmhc()
    if args.renum:
        pdb.clean_docking_count()
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


if __name__ == '__main__':
    main()


