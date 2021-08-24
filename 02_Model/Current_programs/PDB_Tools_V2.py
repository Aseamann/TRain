from math import sqrt
from IGMT_Dict import IGMTDict
import getopt
import sys
import os
import difflib
import shutil
from Bio import PDB
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

class PdbTools:
    # def main(self):
    #     try:
    #         opts, args = getopt.getopt(sys.argv[1:], "i:o:s")
    #     except getopt.GetoptError as err:
    #         sys.stdout = sys.stderrprint(str(err))
    #         usage()
    #         sys.exit(2)
    #     for (opt, arg) in opts:
    #         print("option:", opt)
    #         print("argument:", arg)
    #         if (opt == '-i'):
    #             print(arg)
    #         if (opt == '-o'):
    #             print(arg)
    #         if (opt == '-s'):
    #             print(arg)
    #
    #     # initialize PdbTools
    #     # Takes in a PDB or TXT file
    #     def __init__(self):
    #         self.file_name = sys.argv[1]
    #         self.record_type = ''
    #         self.igmt = IGMTDict('IGMT_Data.txt').get_dict()
    #     # Returns the pdbID of file
    #     def get_pdb_id(self):
    #         with open(self.file_name, 'r') as file:
    #             output = file.readline()
    #             return output[62:66].lower()
    #
    #     print(self.get_pdb_id())
    #
    # if __name__ == "__main__":
    #     main()

    # initialize PdbTools
    def __init__(self, file):
        self.file_name = file
        self.record_type = ''
        self.igmt = IGMTDict('IGMT_Data.txt').get_dict()

    # initialize PdbTools
    def __init__(self):
        self.file_name = ''
        self.record_type = ''
        self.igmt = IGMTDict('IGMT_Data.txt').get_dict()
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

    # Sets record type based on short notation of record
    def set_record_type(self, record_type_in):
        record_type_list = [
            'HEADER', 'OBSLTE', 'TITLE', 'CAVEAT', 'COMPND', 'SOURCE', 'KEYWDS', 'EXPDTA', 'AUTHOR', 'REVDAT',
            'SPRSDE', 'JRNL', 'REMARK', 'DBREF', 'SEQADV', 'SEQRES', 'MODRES', 'HET', 'HETNAME', 'HETSYN',
            'FORMUL', 'HELIX', 'SHEET', 'TURN', 'SSBOND', 'LINK', 'HYDBND', 'SLTBRG', 'CISPEP', 'SITE', 'CRYST1',
            'ORIGX1', 'ORIGX2', 'ORIGX3', 'SCALE1', 'SCALE2', 'SCALE3', 'MTRIX1', 'MTRIX2', 'MTRIX3', 'TVECT',
            'MODEL', 'ATOM', 'SIGATM', 'ANISOU', 'SIGUIJ', 'TER', 'HETATM', 'ENDMDL', 'CONECT', 'MASTER'
        ]
        for i in record_type_list:
            if i == record_type_in:
                self.record_type = record_type_in
                return True
        return False

    # Returns record type in short notation for use in printing strings
    def get_record_type(self):
        return self.record_type

    # Returns record type with 6 char values for searching through PDB file
    def get_record_choice(self):
        record_choice = self.get_record_type()
        length = len(record_choice)
        while length is not 6:
            record_choice = record_choice + ' '
            length += 1
        return record_choice

    # Returns all lines from PDB file with the specified record type
    def record_report(self):
        output = ''
        counter = 0
        flag = False

        # Setting record choice to be used for comparision to file
        record_choice = self.get_record_choice()

        # Reading through file and printing lines with record chosen
        with open(self.file_name, 'r') as file:
            # Prints off first line if the user selected HEADER, couldn't get it to work with for loop below
            if record_choice is 'HEADER':
                file.read(10)
                output += file.readline()
                flag = True
            # Prints off every instance of a line with selected Record Type
            for line in file:
                if record_choice == 'ATOM  ' and line[0:6] == 'TER   ':
                    output += line[6:80]
                    output += '\n'
                    flag = True
                if record_choice == 'ATOM  ' and line[0:6] == 'ATOM  ':
                    output += line[6:80]
                    output += '\n'
                    flag = True
                elif line[0:6] == record_choice:
                    output += line[7:80]
                    output += '\n'
                    flag = True
            # Prints off message to user if file doesn't contain Record Type selected
            if not flag:
                print('File does not contain ' + self.get_record_choice() + ' record type.')
        return output

    # Returns the resolution as an double as reported in REMARK
    def get_resolution(self):
        value = ''
        output = -1
        self.set_record_type('REMARK')
        remark = self.get_record_choice()
        with open(self.file_name, 'r') as file:
            flag = False
            for line in file:
                temp1 = line[0:6]
                if temp1 == remark and not flag:
                    temp = line[6:10]
                    if temp == '   2':
                        flag = True
                elif line[0:6] == remark and flag:
                    value += line[26:29]
                    break
        output = float(value)
        return output

    # Returns the title as a single line
    def get_title(self):
        output = ''
        self.set_record_type('TITLE')
        title = self.get_record_choice()
        with open(self.file_name, 'r') as file:
            for line in file:
                start_pos = 9
                if line[0:6] == title:
                    flag = False
                    temp = line[start_pos:]
                    for char in temp:
                        if char == ' ' and not flag:
                            output += char
                            flag = True
                        elif char == ' ' and flag:
                            break
                        else:
                            output += char
                            flag = False
                start_pos += 1
        return output[1:-1]

    # Returns a list of all chains in PDB file
    def get_chains(self):
        chains = []
        with open(self.file_name, 'r') as file:
            for line in file:
                if line[0:6] == 'ATOM  ':
                    if not chains.__contains__(line[21]):
                        chains.append(line[21])
        return chains

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

    # Finds identical sequence that is contained within a specific chain
    # Need to return specific substring location (not complete)
    def find_seq_in_chain(self, seq_in, chain_in):
        chain = self.get_amino_acid_on_chain(chain_in)
        if chain.__contains__(seq_in):
            return chain.find(seq_in)

    # Returns a dictionary with elements related to specific atom number entered
    def get_atom(self, atom_num):
        with open(self.file_name, 'r') as file:
            for line in file:
                if line[0:6] == 'ATOM  ':
                    if int(line[6:11]) == atom_num and len(line) >= 76:
                        atom = {'atom_num': int(line[6:11]), 'atom_id': line[13:16].strip(), 'atom_comp_id': line[17:20],
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

    # Returns alpha and beta chain IDs based on seq. alignment to 1a07 PDB entry chains. Confirms that it is a
    # partnering TCR chain
    def get_tcr_chains(self):
        matrix = matlist.blosum62
        result = {}
        # Hard coded peptide chains for alpha and beta elements of the TCR_file
        alpha_chain = [
            'KEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFTYRQYSGKSPELIMSIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSDSATYLCAVTTDSTGKLQFGAGTQVVVTPDIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVATSNKSDFACANAFNNSIIPEDTFFPSPESS',
            'QKVTQTQTSISVMEKTTVTMDCVYETQDSSYFLFTYKQTASGEIVFLIRQDSYKKENATVGHYSLNFQKPKSSIGLIITATQIEDSAVYFCAMRGDYGGSGNKLIFGTGTLLSVKP']
        beta_chain = [
            'NAGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSTYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVTEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSTTVNGKEVHSGVSTDPQPLKEQPALNDSRYALSSRLRVSATFTQNPRNHFRCQVQFYGLSENDETTQDRAKPVTQIVSAEATGRAD',
            'VTLLEQNPRTRLVPRGQAVNLRCILKNSQYPTMSTYQQDLQKQLQTLFTLRSPGDKEVKSLPGADYLATRVTDTELRLQVANMSQGRTLYCTCSADRVGNTLYFGEGSRLIV']
        chains = self.get_chains()
        tmp_alpha = []
        tmp_beta = []
        for chain in chains:
            score_alpha = pairwise2.align.globaldx(self.get_amino_acid_on_chain(chain), alpha_chain[0], matrix, score_only=True, penalize_end_gaps=(False, False))
            tmp_alpha.append([float(score_alpha), chain])
            score_beta = pairwise2.align.globaldx(self.get_amino_acid_on_chain(chain), beta_chain[0], matrix, score_only=True, penalize_end_gaps=(False, False))
            tmp_beta.append([float(score_beta), chain])
        alpha = sorted(tmp_alpha)
        beta = sorted(tmp_beta)
        result['ALPHA'] = alpha[-1][1]
        result['BETA'] = beta[-1][1]
        atom = self.first_atom_on_chain(result['ALPHA'])
        position = -1
        while True:
            atom2 = self.first_atom_on_chain(result['BETA'])
            distance1 = self.euclidean_of_atoms(atom['atom_num'], atom2['atom_num'])  # Distance between 1st atom in each chain
            distance2 = self.euclidean_of_atoms(atom['atom_num'], atom2['atom_num'] + 125)  # Distance between 1st and 125 atoms in each chain
            if distance1 <= 46 and abs(distance1 - distance2) <= 15:  # Determines if TCR chains are within typical distance.
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
    def get_tcr_amino_seq_V1(self, tcr_type_in):
        tcr_dict = self.get_tcr_chains()
        for key in tcr_dict:
            if key == tcr_type_in:
                return self.get_amino_acid_on_chain(tcr_dict[key])

    # Checks for gene family of TCR based on IGMT gene family entry of CDR1 and CDR2 sequence
    def get_gene_family_in_tcr(self, chain):
        output = {}
        tcr_seq = self.get_tcr_amino_seq_V1(chain)
        for gene_family in self.igmt:
            if tcr_seq.__contains__(self.igmt[gene_family]['CDR1-IMGT']):
                start = tcr_seq.find(self.igmt[gene_family]['CDR1-IMGT'])
                if tcr_seq[start:].__contains__(self.igmt[gene_family]['CDR2-IMGT']):
                    output[gene_family] = self.igmt[gene_family]
        return output

    # Previous version, saved for any issues arising in new method.
    def get_gene_family_in_tcr_save(self, chain):
        output = {}
        tcr_seq = self.get_tcr_amino_seq_V1(chain)
        for gene_family in self.igmt:
            for cdr_type in self.igmt[gene_family]:
                if cdr_type.__contains__('CDR1') or cdr_type.__contains__('CDR2'):
                    if tcr_seq.__contains__(self.igmt[gene_family][cdr_type]):
                        print(self.igmt[gene_family][cdr_type])
                        if len(self.igmt[gene_family][cdr_type]) > 3:
                            output[gene_family] = self.igmt[gene_family]
        return output

    # Returns an array of atoms located on specific cdr chain based on cdr type and alpha or beta chain
    # Chain is either 'ALPHA' or 'BETA'
    def get_atoms_on_cdr(self, cdr_type_in, chain):
        atoms = []
        cdr_type = ''
        if cdr_type_in == 1:
            cdr_type = 'CDR1-IMGT'
        elif cdr_type_in == 2:
            cdr_type = 'CDR2-IMGT'
        tcr_chain = self.get_tcr_amino_seq_V1(chain)
        tcr_dict = self.get_gene_family_in_tcr(chain)
        print(tcr_dict)
        dict_key_1 = next(iter(tcr_dict))
        comp_id_start = tcr_chain.find(tcr_dict[dict_key_1][cdr_type])
        comp_id_end = comp_id_start + len(tcr_dict[dict_key_1][cdr_type])
        self.set_record_type('ATOM')
        file = self.record_report().splitlines()
        flag = True
        for line in file:
            if line[15] == self.get_tcr_chains()[chain]:
                number = int(line[16:20])
                if flag:
                    comp_id_start += number
                    comp_id_end += number - 1
                    flag = False
                if number <= comp_id_end:
                    if comp_id_start < number:
                        comp_id_start += 1
                    if comp_id_start == number:
                        atoms.append(self.get_atom(int(line[0:5])))
                else:
                    break
        print(atoms)
        return atoms

    def get_cdr_seq(self, cdr_type_in, chain):
        cdr_type = ''
        if cdr_type_in == 1:
            cdr_type = 'CDR1-IMGT'
        elif cdr_type_in == 2:
            cdr_type = 'CDR2-IMGT'
        tcr_chain = self.get_tcr_amino_seq_V1(chain)
        print(tcr_chain)
        tcr_dict = self.get_gene_family_in_tcr(chain)
        print(tcr_dict)
        if bool(tcr_dict):
            dict_key_1 = next(iter(tcr_dict))
            comp_id_start = tcr_chain.find(tcr_dict[dict_key_1][cdr_type])
            comp_id_end = comp_id_start + len(tcr_dict[dict_key_1][cdr_type])
            return tcr_dict[dict_key_1][cdr_type]

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
        self.set_record_type('SEQRES')
        file_seqres = self.record_report().splitlines()
        lowest = 1000
        chain = ''
        for line in file_seqres:
            if int(line[7:10]) < lowest:
                lowest = int(line[7:10])
                chain = line[4]
        return chain

    # Returns the MHC chain - based on COMPND section labeling HISTOCOMPATIBILITY ANTIGEN
    # WILL LATER CONVERT THIS TO HOW I CHECK FOR TCR CHAINS
    def get_mhc_chain(self):
        self.set_record_type('COMPND')
        file = self.record_report().splitlines()
        chains = self.get_chains()
        flag = False
        for line in file:
            if flag and line[11] in chains:
                return line[11]
            elif line[4:12] == 'MOLECULE' and line.__contains__('HISTOCOMPATIBILITY ANTIGEN'):
                flag = True

    # Returns the atoms on the MHC chain = based on COMPND section labeling HISTOCOMPATIBILITY ANTIGEN
    def get_atoms_on_mhc(self):
        atoms = []
        self.set_record_type('ATOM')
        file_atoms = self.record_report().splitlines()
        chain = self.get_mhc_chain()
        for line in file_atoms:
            if line[15] == chain:
                atoms.append(self.get_atom(int(line[0:5])))
        atoms.pop(len(atoms) - 1)
        return atoms

    # Returns atoms in TCR chain selected that are within 4.5 angstroms of peptide, based on either alpha or beta
    def get_atoms_tcr_to_peptide(self, chain):
        tcr = []
        tcr_2 = []
        cdr_1_atoms = self.get_atoms_on_cdr(1, chain)
        cdr_2_atoms = self.get_atoms_on_cdr(2, chain)
        peptide_atoms = self.get_atoms_on_peptide()
        for atom_p in peptide_atoms:
            for atom_1 in cdr_1_atoms:
                if self.euclidean_of_atoms(atom_1['atom_num'], atom_p['atom_num']) <= 4.5:
                    tcr.append([1, atom_1['atom_num'], atom_p['atom_num']])
                    break
            for atom_2 in cdr_2_atoms:
                if self.euclidean_of_atoms(atom_2['atom_num'], atom_p['atom_num']) <= 4.5:
                    tcr_2.append([2, atom_2['atom_num'], atom_p['atom_num']])
                    break
        tcr.extend(tcr_2)
        return tcr

    # Returns atoms within the TCR chain selected that are within 4.5 angstroms of MHC, based on either beta or alpha
    def get_atoms_tcr_to_mhc(self, chain):
        tcr = []
        tcr_2 = []
        cdr_1_atoms = self.get_atoms_on_cdr(1, chain)
        cdr_2_atoms = self.get_atoms_on_cdr(2, chain)
        mhc_atoms = self.get_atoms_on_mhc()
        count = 180
        for atom_mhc in mhc_atoms:
            for atom_1 in cdr_1_atoms:
                if self.euclidean_of_atoms(atom_1['atom_num'], atom_mhc['atom_num']) <= 4.5:
                    tcr.append([1, atom_1['atom_num'], atom_mhc['atom_num']])
            for atom_2 in cdr_2_atoms:
                if self.euclidean_of_atoms(atom_2['atom_num'], atom_mhc['atom_num']) <= 4.5:
                    tcr_2.append([2, atom_2['atom_num'], atom_mhc['atom_num']])
            count += 1
            if atom_mhc['comp_num'] >= count:
                break
        tcr.extend(tcr_2)
        return tcr

    # Returns plain text of TCR chains with PBD/chain for partnering TCR chains selected.
    def get_pisces_entry_tcr(self, logfile):
        tcr_list = self.get_tcr_chains()
        with open(logfile, 'a+') as f:
            for key in tcr_list:
                line = self.get_pdb_id() + tcr_list[key] + '\n'
                f.write(line)

    # Returns plain text of ALL chains with PBD/chain.
    def get_pisces_entry_all(self, piscesfile):
        chain_list = self.get_chains()
        with open(piscesfile, 'a+') as f:
            for each in chain_list:
                line = str(self.get_pdb_id()) + str(each) + '\n'
                f.write(str(line))

    # Returns a reformatted PDB with just the TCR atom cord. and has relabeled chains with ALPHA = A; BETA = B
    def clean_tcr(self, dir_start='****'):
        if dir_start != '****':
            tcr = dir_start + '%s_tcr.pdb' % (self.get_pdb_id())
        else:
            tcr = '%s_tcr.pdb' % (self.get_pdb_id())
        tcr_list = self.get_tcr_chains()
        atom_count = 0
        flag = False
        with open(self.file_name) as f:
            with open(tcr, 'w+') as f1:
                for line in f:
                    if line[0:6] == 'HEADER':
                        f1.write(line)
                        flag = True
                    if flag:
                        f1.write('EXPDTA    THEORETICAL MODEL    CLEAN TCR ALPHA:A BETA:B\n')
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
                                f1.write(line.replace(num, str(atom_count).rjust(5), 1))

    def recount_tcr(self):
        atom_count = 0
        flag_a = False
        flag_b = False
        res_alpha_count = 1
        res_beta_count = 1
        previous_count_a = -1
        previous_count_b = -1
        with open(self.file_name, 'r') as f:
            tcr = f.readlines()
        with open(self.file_name, 'w+') as f1:
            for line in tcr:
                if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                    if line[16] != 'B' and line[26] == ' ':
                        if line[21] == "A" or line[21] == "B":
                            if line[21] == "A":
                                if int(line[22:26]) != previous_count_a and not flag_a:
                                    previous_count_a = int(line[22:26])
                                    flag_a = True
                                elif int(line[22:26]) != previous_count_a and flag_a:
                                    previous_count_a = int(line[22:26])
                                    res_alpha_count += 1
                                line = line[:21] + 'A' + line[22:]
                                line = line[:22] + str(res_alpha_count).rjust(4) + line[26:]
                            elif line[21] == "B":
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
                            f1.write(line.replace(num, str(atom_count).rjust(5), 1))
                else:
                    f1.write(line)

    # Returns a reformatted PDB with TCRpMHC complex renumbered for RosettaDock
    # New count on amino acids so each chain starts at previous count
    def clean_docking_count(self, rename='****'):
        if rename != '****':
            tcr = rename
        else:
            tcr = 'Docking/test_docking.pdb'
        tcr_list = self.get_tcr_chains()
        print(tcr_list)
        atom_count = 0
        flag = False
        flag_a = False
        flag_b = False
        flag_x = False
        flag_res = False
        last_res = ''
        res_alpha_count = 1
        alpha_cut = 107
        res_beta_count = 1
        beta_cut = 113
        previous_count_a = -1
        previous_count_b = -1
        previous_count_x = -1
        res_count = 1
        with open(self.file_name) as f:
            with open(tcr, 'w+') as f1:
                for line in f:
                    if line[0:6] == 'HEADER':
                        f1.write(line)
                        flag = True
                    if flag:
                        f1.write('EXPDTA    THEORETICAL MODEL    CLEAN TCR ALPHA:A BETA:B\n')
                        flag = False
                    if line[0:6] == 'ATOM  ':  # Only write over atoms
                        if line[16] != 'B' and line[26] == ' ':  # Don't allow secondary atoms
                            if last_res != line[21] and flag_res:
                                f1.write('TER\n')
                            elif last_res != line[21] and not flag_res:
                                flag_res = True
                            if line[21] == tcr_list.get('ALPHA') and res_alpha_count <= alpha_cut:
                                if int(line[22:26]) != previous_count_a and not flag_a:
                                    previous_count_a = int(line[22:26])
                                    flag_a = True
                                elif int(line[22:26]) != previous_count_a and flag_a:
                                    previous_count_a = int(line[22:26])
                                    res_alpha_count += 1
                                    res_count += 1
                                line = line[:22] + str(res_count).rjust(4) + line[26:]
                                num = line[6:11]
                                atom_count += 1
                                last_res = line[21]
                                if line[16] == 'A':
                                    line = line[:16] + ' ' + line[17:]
                                if res_alpha_count <= alpha_cut:
                                    f1.write(line.replace(num, str(atom_count).rjust(5), 1))
                            elif line[21] == tcr_list.get('BETA') and res_beta_count <= beta_cut:
                                if int(line[22:26]) != previous_count_b and not flag_b:
                                    previous_count_b = int(line[22:26])
                                    flag_b = True
                                elif int(line[22:26]) != previous_count_b and flag_b:
                                    previous_count_b = int(line[22:26])
                                    res_beta_count += 1
                                    res_count += 1
                                line = line[:22] + str(res_count).rjust(4) + line[26:]
                                num = line[6:11]
                                atom_count += 1
                                last_res = line[21]
                                if line[16] == 'A':
                                    line = line[:16] + ' ' + line[17:]
                                if res_beta_count <= beta_cut:
                                    f1.write(line.replace(num, str(atom_count).rjust(5), 1))
                            elif line[21] != tcr_list.get('BETA') and line[21] != tcr_list.get('ALPHA'):
                                if res_count == 179:
                                    print("test")
                                if int(line[22:26]) != previous_count_x and not flag_x:
                                    previous_count_x = int(line[22:26])
                                    flag_x = True
                                elif int(line[22:26]) != previous_count_x and flag_x:
                                    previous_count_x = int(line[22:26])
                                    res_count += 1
                                line = line[:22] + str(res_count).rjust(4) + line[26:]
                                num = line[6:11]
                                atom_count += 1
                                last_res = line[21]
                                if line[16] == 'A':
                                    line = line[:16] + ' ' + line[17:]
                                if res_beta_count <= beta_cut:
                                    f1.write(line.replace(num, str(atom_count).rjust(5), 1))
                f1.write('TER\nEND\n')

    # Returns a reformatted PDB with just the TCR atom cord. and has relabeled chains with ALPHA = A; BETA = B
    # New count on amino acids so each chain starts at 1
    def clean_tcr_count(self, dir_start='****'):
        if dir_start != '****':
            tcr = dir_start + '%s_tcr.pdb' % (self.get_pdb_id())
        else:
            tcr = '%s_tcr.pdb' % (self.get_pdb_id())
        tcr_list = self.get_tcr_chains()
        atom_count = 0
        flag = False
        flag_a = False
        flag_b = False
        res_alpha_count = 1
        res_beta_count = 1
        previous_count_a = -1
        previous_count_b = -1
        res_count = 0
        with open(self.file_name) as f:
            with open(tcr, 'w+') as f1:
                for line in f:
                    if line[0:6] == 'HEADER':
                        f1.write(line)
                        flag = True
                    if flag:
                        f1.write('EXPDTA    THEORETICAL MODEL    CLEAN TCR ALPHA:A BETA:B\n')
                        flag = False
                    if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                        if line[16] != 'B' and line[26] == ' ':
                            if line[21] == tcr_list.get('ALPHA') or line[21] == tcr_list.get('BETA'):
                                if line[21] == tcr_list.get('ALPHA'):
                                    if int(line[22:26]) != previous_count_a and not flag_a:
                                        previous_count_a = int(line[22:26])
                                        flag_a = True
                                    elif int(line[22:26]) != previous_count_a and flag_a:
                                        previous_count_a = int(line[22:26])
                                        res_alpha_count += 1
                                    line = line[:21] + 'A' + line[22:]
                                    line = line[:22] + str(res_alpha_count).rjust(4) + line[26:]
                                elif line[21] == tcr_list.get('BETA'):
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
                                f1.write(line.replace(num, str(atom_count).rjust(5), 1))

    # Returns a reformatted PDB with just the TCR atom cord. and has relabeled chains with ALPHA = A; BETA = B
    # Also uses trimmed TCR seq.
    # New count on amino acids so each chain starts at 1
    def clean_tcr_count_trim(self, dir_start='****'):
        if dir_start != '****':
            tcr = dir_start + '%s_tcr.pdb' % (self.get_pdb_id())
        else:
            tcr = '%s_tcr.pdb' % (self.get_pdb_id())
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
        with open(self.file_name) as f:
            with open(tcr, 'w+') as f1:
                for line in f:
                    if line[0:6] == 'HEADER':
                        f1.write(line)
                        flag = True
                    if flag:
                        f1.write('EXPDTA    THEORETICAL MODEL    CLEAN TCR ALPHA:A BETA:B\n')
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
                                        f1.write(line.replace(num, str(atom_count).rjust(5), 1))
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
                                        f1.write(line.replace(num, str(atom_count).rjust(5), 1))
        print(res_alpha_count)
        print(res_beta_count)

    # Creates a new PDB file containing only information needed for TCR for the TCR of the original PDB file
    def spilt_tcr(self):
        tcr = '%s_tcr.pdb' % (self.get_pdb_id())
        tcr_list = self.get_tcr_chains()
        helix_count = 0
        sheet_count = 0
        atom_count = 0
        compare_conect = {}
        with open(self.file_name) as f:
            with open(tcr, 'w+') as f1:
                for line in f:
                    left_conect = 6
                    right_conect = 11
                    if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                        if line[16] != 'B':
                            if line[21] == tcr_list.get('ALPHA') or line[21] == tcr_list.get('BETA'):
                                num = line[6:11]
                                atom_count += 1
                                f1.write(line.replace(num, str(atom_count).rjust(5), 1))
                                compare_conect[num] = atom_count
                    elif line[0:6] == 'HELIX ':
                        if line[19] == tcr_list.get('ALPHA') or line[19] == tcr_list.get('BETA'):
                            num = line[6:10]
                            helix_count += 1
                            f1.write(line.replace(num, str(helix_count).rjust(4), 1))
                    elif line[0:6] == 'SHEET ':
                        if line[21] == tcr_list.get('ALPHA') or line[21] == tcr_list.get('BETA'):
                            num = line[6:10]
                            sheet_count += 1
                            f1.write(line.replace(num, str(sheet_count).rjust(4), 1))
                    # elif line[0:6] == 'HETATM':
                        # if line[21] == tcr_list.get('ALPHA') or line[21] == tcr_list.get('BETA'):
                        #     num = line[6:11]
                        #     atom_count += 1
                        #     f1.write(line.replace(num, str(atom_count).rjust(5), 1))
                    elif line[0:6] == 'CONECT':
                        if compare_conect.__contains__(line[left_conect:right_conect]):
                            while compare_conect.__contains__(line[left_conect:right_conect]):
                                line_update = line.replace(line[left_conect:right_conect]
                                        ,str(compare_conect[line[left_conect:right_conect]]).rjust(5), 1)
                                left_conect += 5
                                right_conect += 5
                            f1.write(line_update)
                    else:
                        if line[0:6] != 'MASTER' and line[0:6] != 'HETATM' and line[0:6] != 'ANISOU':
                            f1.write(line)

    # Creates a new PDB file with information for only the peptide of the original PDB file
    def spilt_p(self):
        tcr = '%s_peptide.pdb' % (self.get_pdb_id())
        peptide = self.get_peptide_chain()
        helix_count = 0
        sheet_count = 0
        atom_count = 0
        compare_conect = {}
        with open(self.file_name) as f:
            with open(tcr, 'w+') as f1:
                for line in f:
                    left_conect = 6
                    right_conect = 11
                    if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                        if line[21] == peptide:
                            num = line[6:11]
                            atom_count += 1
                            f1.write(line.replace(num, str(atom_count).rjust(5), 1))
                            compare_conect[num] = atom_count
                    elif line[0:6] == 'HELIX ':
                        if line[19] == peptide:
                            num = line[6:10]
                            helix_count += 1
                            f1.write(line.replace(num, str(helix_count).rjust(4), 1))
                    elif line[0:6] == 'SHEET ':
                        if line[21] == peptide:
                            num = line[6:10]
                            sheet_count += 1
                            f1.write(line.replace(num, str(sheet_count).rjust(4), 1))
                    elif line[0:6] == 'HETATM':
                        if line[21] == peptide:
                            num = line[6:11]
                            atom_count += 1
                            f1.write(line.replace(num, str(atom_count).rjust(5), 1))
                    elif line[0:6] == 'CONECT':
                        if compare_conect.__contains__(line[left_conect:right_conect]):
                            while compare_conect.__contains__(line[left_conect:right_conect]):
                                line_update = line.replace(line[left_conect:right_conect]
                                                           ,
                                                           str(compare_conect[line[left_conect:right_conect]]).rjust(5),
                                                           1)
                                left_conect += 5
                                right_conect += 5
                            f1.write(line_update)
                    else:
                        if line[0:6] != 'MASTER':
                            f1.write(line)
    # Creates a new PDB file with information for only the MHC of the original PDB file
    def spilt_mhc(self):
        tcr = '%s_mhc.pdb' % (self.get_pdb_id())
        mhc = self.get_mhc_chain()
        helix_count = 0
        sheet_count = 0
        atom_count = 0
        compare_conect = {}
        with open(self.file_name) as f:
            with open(tcr, 'w+') as f1:
                for line in f:
                    left_conect = 6
                    right_conect = 11
                    if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                        if line[21] == mhc:
                            num = line[6:11]
                            atom_count += 1
                            f1.write(line.replace(num, str(atom_count).rjust(5), 1))
                            compare_conect[num] = atom_count
                    elif line[0:6] == 'HELIX ':
                        if line[19] == mhc:
                            num = line[6:10]
                            helix_count += 1
                            f1.write(line.replace(num, str(helix_count).rjust(4), 1))
                    elif line[0:6] == 'SHEET ':
                        if line[21] == mhc:
                            num = line[6:10]
                            sheet_count += 1
                            f1.write(line.replace(num, str(sheet_count).rjust(4), 1))
                    elif line[0:6] == 'HETATM':
                        if line[21] == mhc:
                            num = line[6:11]
                            atom_count += 1
                            f1.write(line.replace(num, str(atom_count).rjust(5), 1))
                    elif line[0:6] == 'CONECT':
                        if compare_conect.__contains__(line[left_conect:right_conect]):
                            while compare_conect.__contains__(line[left_conect:right_conect]):
                                line_update = line.replace(line[left_conect:right_conect]
                                                           ,
                                                           str(compare_conect[line[left_conect:right_conect]]).rjust(5),
                                                           1)
                                left_conect += 5
                                right_conect += 5
                            f1.write(line_update)
                    else:
                        if line[0:6] != 'MASTER':
                            f1.write(line)

    # Creates a new PDB file containing only information needed for TCR for the TCR of the original PDB file
    def spilt_tcr_alpha(self):
        tcr = '%s_tcr.pdb' % (self.get_pdb_id())
        tcr_list = self.get_tcr_chains()
        helix_count = 0
        sheet_count = 0
        atom_count = 0
        compare_conect = {}
        with open(self.file_name) as f:
            with open(tcr, 'w+') as f1:
                for line in f:
                    left_conect = 6
                    right_conect = 11
                    if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                        if line[16] != 'B':
                            if line[21] == tcr_list.get('ALPHA'):
                                num = line[6:11]
                                atom_count += 1
                                f1.write(line.replace(num, str(atom_count).rjust(5), 1))
                                compare_conect[num] = atom_count
                    elif line[0:6] == 'HELIX ':
                        if line[19] == tcr_list.get('ALPHA'):
                            num = line[6:10]
                            helix_count += 1
                            f1.write(line.replace(num, str(helix_count).rjust(4), 1))
                    elif line[0:6] == 'SHEET ':
                        if line[21] == tcr_list.get('ALPHA') :
                            num = line[6:10]
                            sheet_count += 1
                            f1.write(line.replace(num, str(sheet_count).rjust(4), 1))
                    # elif line[0:6] == 'HETATM':
                    # if line[21] == tcr_list.get('ALPHA') or line[21] == tcr_list.get('BETA'):
                    #     num = line[6:11]
                    #     atom_count += 1
                    #     f1.write(line.replace(num, str(atom_count).rjust(5), 1))
                    elif line[0:6] == 'CONECT':
                        if compare_conect.__contains__(line[left_conect:right_conect]):
                            while compare_conect.__contains__(line[left_conect:right_conect]):
                                line_update = line.replace(line[left_conect:right_conect]
                                                           ,
                                                           str(compare_conect[line[left_conect:right_conect]]).rjust(5),
                                                           1)
                                left_conect += 5
                                right_conect += 5
                            f1.write(line_update)
                    else:
                        if line[0:6] != 'MASTER' and line[0:6] != 'HETATM' and line[0:6] != 'ANISOU':
                            f1.write(line)

    # Returns a formatted pir file for a specific PDB for the ALPHA chain
    # Removed _alpha from writes
    def pir_alpha(self):
        # pir_name = '%s_alpha.pir' % (self.get_pdb_id())
        pir_name = 'ALL_ALPHA_PIR.pir'
        tcr_list = self.get_tcr_chains()
        tcr_alpha_chain = self.get_tcr_amino_seq_V1('ALPHA')
        chain_length = len(tcr_alpha_chain)
        count = 1
        with open(pir_name, 'a+') as f:
            f.write('C; Produced by PDB_Tools_V2\n\n')
            f.write('>P1;' + self.get_pdb_id() + tcr_list['ALPHA'] + '\n')
            f.write('structureX:' + self.get_pdb_id() + ':   1 :' + tcr_list['ALPHA'] + ': ' + str(len(tcr_alpha_chain))
                    + ' :' + tcr_list['ALPHA']
                    + ': : :' + str(self.get_resolution()).rjust(5) + ":-1.00\n")
            for aa in tcr_alpha_chain:
                if count % 76 != 0:
                    f.write(aa)
                    count += 1
                else:
                    f.write('\n' + aa)
                    count += 1
            f.write('*\n')

    # Returns a formatted pir file for a specific PDB for the BETA chain
    def pir_beta(self):
        # pir_name = '%s_beta.pir' % (self.get_pdb_id())
        pir_name = 'ALL_BETA_PIR.pir'
        tcr_list = self.get_tcr_chains()
        tcr_beta_chain = self.get_tcr_amino_seq_V1('BETA')
        count = 1
        with open(pir_name, 'a+') as f:
            f.write('>P1;' + self.get_pdb_id() + '_beta\n')
            f.write('structureX:' + self.get_pdb_id() + ':FIRST:' + tcr_list['BETA'] + ':LAST:' + tcr_list['BETA']
                    + ': : :' + str(self.get_resolution()) + ":-1.00\n")
            for aa in tcr_beta_chain:
                if count % 76 != 0:
                    f.write(aa)
                    count += 1
                else:
                    f.write('\n' + aa)
                    count += 1
            f.write('*\n')

    def pir_tcr(self):
        pir_name = 'ALL_TCR_PIR.pir'
        tcr_list = self.get_tcr_chains()
        tcr_alpha_chain = self.get_tcr_amino_seq_V1('ALPHA')
        tcr_beta_chain = self.get_tcr_amino_seq_V1('BETA')
        tcr_combined = tcr_beta_chain + '/' + tcr_alpha_chain
        count = 1
        with open(pir_name, 'a+') as f:
            f.write('C; Produced by PDB_Tools_V2\n\n')
            f.write('>P1;' + self.get_pdb_id() + '\n')
            f.write('structureX:' + self.get_pdb_id() + ':FIRST:' + tcr_list['BETA'] + ':LAST:' + tcr_list['ALPHA']
                    + ': : :' + str(self.get_resolution()) + ":-1.00\n")
            for aa in tcr_combined:
                if count % 76 != 0:
                    f.write(aa)
                    count += 1
                else:
                    f.write('\n' + aa)
                    count += 1
            f.write('*\n')

    # Returns a formatted ali file for a specific PDB for the ALPHA chain
    def ali_alpha(self):
        ali_name = '%s_alpha.ali' % (self.get_pdb_id())
        tcr_list = self.get_tcr_chains()
        tcr_alpha_chain = self.get_tcr_amino_seq_V1('ALPHA')[0:115]
        chain_length = len(tcr_alpha_chain)
        count = 1
        with open(ali_name, 'a+') as f:
            f.write('C; Produced by PDB_Tools_V2\n\n')
            f.write('>P1;' + self.get_pdb_id() + tcr_list['ALPHA'] + '\n')
            f.write('sequence:' + self.get_pdb_id() + tcr_list['ALPHA'] + ':::::::0.00: 0.00\n')
            for aa in tcr_alpha_chain:
                if count % 76 != 0:
                    f.write(aa)
                    count += 1
                else:
                    f.write('\n' + aa)
                    count += 1
            f.write('*\n')

    # Runs ali_alpha but returns as a string rather than a new file as ali_alpha does.
    def ali_alpha_text(self):
        output = ''
        tcr_list = self.get_tcr_chains()
        tcr_alpha_chain = self.get_tcr_amino_seq_V1('ALPHA')[0:115]
        count = 1
        output += 'C; Produced by PDB_Tools_V2\n\n'
        output += '>P1;' + self.get_pdb_id() + tcr_list['ALPHA'] + '\n'
        output += 'sequence:' + self.get_pdb_id() + tcr_list['ALPHA'] + ':::::::0.00: 0.00\n'
        for aa in tcr_alpha_chain:
            if count % 76 != 0:
                output += aa
                count += 1
            else:
                output += '\n' + aa
                count += 1
        output += '*\n'
        return output

    # Returns a formatted ali file for a specific PDB for the BETA chain
    def ali_beta(self):
        ali_name = '%s_beta.ali' % (self.get_pdb_id())
        tcr_list = self.get_tcr_chains()[0:115]
        tcr_beta_chain = self.get_tcr_amino_seq_V1('BETA')
        count = 1
        with open(ali_name, 'a+') as f:
            f.write('C; Produced by PDB_Tools_V2\n\n')
            f.write('>P1;' + self.get_pdb_id() + tcr_list['BETA'] + '\n')
            f.write('sequence:' + self.get_pdb_id() + tcr_list['BETA'] + ':::::::0.00: 0.00\n')
            for aa in tcr_beta_chain:
                if count % 76 != 0:
                    f.write(aa)
                    count += 1
                else:
                    f.write('\n' + aa)
                    count += 1
            f.write('*\n')

    # Returns a formatted ali file for a specific PDB for both Alpha and Beta chain
    def ali_tcr(self):
        ali_name = '%s_tcr.ali' % (self.get_pdb_id())
        tcr_list = self.get_tcr_chains()
        tcr_alpha_chain = self.get_tcr_amino_seq_V1('ALPHA')
        tcr_beta_chain = self.get_tcr_amino_seq_V1('BETA')
        tcr_combined = tcr_beta_chain + '/' + tcr_alpha_chain
        count = 1
        with open(ali_name, 'a+') as f:
            f.write('C; Produced by PDB_Tools_V2\n\n')
            f.write('>P1;' + self.get_pdb_id() + '\n')
            f.write('sequence:' + self.get_pdb_id() + ':::::::0.00: 0.00\n')
            for aa in tcr_combined:
                if count % 76 != 0:
                    f.write(aa)
                    count += 1
                else:
                    f.write('\n' + aa)
                    count += 1
            f.write('*\n')

    # Returns a fasta formatted file of PDB file submitted. Default to both chains or input for specific
    def fasta_TCR(self, file_name='result.fasta', chain='*****'):
        tcr_alpha_chain = self.get_amino_acid_on_chain('A')
        tcr_beta_chain = self.get_amino_acid_on_chain('B')
        pdb_id = self.get_pdb_id()
        count_1 = 1
        count_2 = 1
        count_3 = 1
        with open(file_name, 'a+') as f:
            if chain == 'ALPHA':
                if tcr_alpha_chain != '':
                    f.write('>' + pdb_id + 'A\n')
                    for aa in tcr_alpha_chain:
                        if count_1 % 81 != 0:
                            f.write(aa)
                            count_1 += 1
                        else:
                            f.write('\n' + aa)
                            count_1 = 2
                    f.write('\n')
            elif chain == 'BETA':
                if tcr_beta_chain != '':
                    f.write('>' + pdb_id + 'B\n')
                    for aa in tcr_beta_chain:
                        if count_2 % 81 != 0:
                            f.write(aa)
                            count_2 += 1
                        else:
                            f.write('\n' + aa)
                            count_2 = 2
                    f.write('\n')
            elif chain == '*****':
                if tcr_alpha_chain != '' and tcr_beta_chain != '':
                    f.write('>' + pdb_id + '\n')
                    tcr = tcr_alpha_chain + tcr_beta_chain
                    for aa in tcr:
                        if count_3 % 81 != 0:
                            f.write(aa)
                            count_3 += 1
                        else:
                            f.write('\n' + aa)
                            count_3 = 2
                    f.write('\n')

    # Returns a formatted pir file for a specific PDB for the BETA chain
    def pir_tcr_single(self):
        pir_name = '%s_tcr.pir' % (self.get_pdb_id())
        tcr_list = self.get_tcr_chains()
        tcr_alpha_chain = self.get_tcr_amino_seq_V1('ALPHA')
        tcr_beta_chain = self.get_tcr_amino_seq_V1('BETA')
        tcr_combined = tcr_beta_chain + '/' + tcr_alpha_chain
        count = 1
        with open(pir_name, 'a+') as f:
            f.write('>P1;' + self.get_pdb_id() + '\n')
            f.write('structureX:' + self.get_pdb_id() + ':FIRST:' + tcr_list['BETA'] + ':LAST:' + tcr_list['ALPHA']
                    + ': : :' + str(self.get_resolution()) + ":-1.00\n")
            for aa in tcr_combined:
                if count % 76 != 0:
                    f.write(aa)
                    count += 1
                else:
                    f.write('\n' + aa)
                    count += 1
            f.write('*\n')

    # Groups PDB files based on gene family from IGMT Dict, paths are setup for use on my computer only
    def group_gene_family(self):
        try:
            family_name_alpha = self.get_gene_family_in_tcr('ALPHA')
            print(family_name_alpha)
            family_name_beta = self.get_gene_family_in_tcr('BETA')
            print(family_name_beta)
            first_gene_alpha = iter(family_name_alpha.keys())
            first_gene_beta = iter(family_name_beta.keys())
            newpath_alpha = r'/Users/austinseamann/PycharmProjects/TCRDOCK/GENE_FAMILY/%s' % next(first_gene_alpha)
            newpath_beta = r'/Users/austinseamann/PycharmProjects/TCRDOCK/GENE_FAMILY/%s' % next(first_gene_beta)
            newfile_alpha = newpath_alpha + '/%s_alpha.pdb' % self.get_pdb_id()
            newfile_beta = newpath_beta + '/%s_beta.pdb' % self.get_pdb_id()
            if not os.path.exists(newpath_alpha):
                os.makedirs(newpath_alpha)
            if not os.path.exists(newpath_beta):
                os.makedirs(newpath_beta)
            with open(self.file_name, 'r') as f:
                with open(newfile_alpha, 'w+') as f1:
                    for line in f:
                        f1.write(line)
            with open(self.file_name, 'r') as f2:
                with open(newfile_beta, 'w+') as f3:
                    for line1 in f2:
                        f3.write(line1)
        except:
            print('No Gene Family Found')

    def rmsd_all(self, sample):
        rmsd = {}
        rmsd['ALL'] = self.rmsd(sample, 'ALL', 'ALL')
        rmsd['ALL_CDR1'] = self.rmsd(sample, 'ALL', 'CDR1')
        rmsd['ALL_CDR2'] = self.rmsd(sample, 'ALL', 'CDR2')
        # rmsd['ALL_CDR3'] = self.rmsd(sample, 'ALL', 'CDR3')
        rmsd['ALPHA'] = self.rmsd(sample, 'ALPHA', 'ALL')
        rmsd['ALPHA_CDR1'] = self.rmsd(sample, 'ALPHA', 'CDR1')
        rmsd['ALPHA_CDR2'] = self.rmsd(sample, 'ALPHA', 'CDR2')
        # rmsd['ALPHA_CDR3'] = self.rmsd(sample, 'ALPHA', 'CDR3')
        rmsd['BETA'] = self.rmsd(sample, 'BETA', 'ALL')
        # rmsd['BETA_CDR1'] = self.rmsd(sample, 'BETA', 'CDR1')
        # rmsd['BETA_CDR2'] = self.rmsd(sample, 'BETA', 'CDR2')
        # rmsd['BETA_CDR3'] = self.rmsd(sample, 'BETA', 'CDR3')
        return rmsd

    def rmsd(self, sample, chain='****', par='****'):
        length_alpha = len(self.get_amino_acid_on_chain('A'))
        length_beta = len(self.get_amino_acid_on_chain('B'))
        ref = self.get_file_name()
        self.set_file_name(sample)
        tcr_chain = self.get_tcr_chains()
        if length_alpha > len(self.get_amino_acid_on_chain(tcr_chain['ALPHA'])):
            length_alpha = len(self.get_amino_acid_on_chain(tcr_chain['ALPHA']))
        if length_beta > len(self.get_amino_acid_on_chain(tcr_chain['BETA'])):
            length_beta = len(self.get_amino_acid_on_chain(tcr_chain['BETA']))
        self.set_file_name(ref)
        if chain == 'ALL' or chain == '****':
            if par == 'ALL':
                return self.rmsd_helper(sample, 'ALL', 0, length_alpha, 0, length_beta)
            elif par == 'CDR1':
                return self.rmsd_helper(sample, 'ALL', 27, 38, 27, 38)
            elif par == 'CDR2':
                return self.rmsd_helper(sample, 'ALL', 56, 65, 56, 65)
            elif par == 'CDR3':
                return self.rmsd_helper(sample, 'ALL', 105, 116, 105, 116)
        elif chain == 'ALPHA':
            if par == 'ALL':
                return self.rmsd_helper(sample, 'ALPHA', 0, length_alpha, 0, 0)
            elif par == 'CDR1':
                return self.rmsd_helper(sample, 'ALPHA', 27, 38, 0, 0)
            elif par == 'CDR2':
                return self.rmsd_helper(sample, 'ALPHA', 56, 65, 0, 0)
            elif par == 'CDR3':
                return self.rmsd_helper(sample, 'ALPHA', 105, 116, 0, 0)
        elif chain == 'BETA':
            if par == 'ALL':
                return self.rmsd_helper(sample, 'BETA', 0, 0, 0, length_beta)
            elif par == 'CDR1':
                return self.rmsd_helper(sample, 'BETA', 0, 0, 27, 38)
            elif par == 'CDR2':
                return self.rmsd_helper(sample, 'BETA', 0, 0, 56, 65)
            elif par == 'CDR3':
                return self.rmsd_helper(sample, 'BETA', 0, 0, 105, 116)

    def rmsd_helper(self, sample, chain, start_a=0, end_a=0, start_b=0, end_b=0):
        chains_ref = self.get_tcr_chains()
        temp = self.get_file_name()
        self.set_file_name(sample)
        chains_sample = self.get_tcr_chains()
        self.set_file_name(temp)
        if chain == 'ALL':
            chain_id_ref = 'ALL'
            chain_id_sample = 'ALL'
        else:
            chain_id_ref = chains_ref[chain]
            chain_id_sample = chains_sample[chain]
        pdb_parser = PDB.PDBParser(QUIET=True)

        ref_structure = pdb_parser.get_structure("reference", self.get_file_name())
        sample_structure = pdb_parser.get_structure("sample", sample)

        ref_model = ref_structure[0]
        sample_model = sample_structure[0]

        ref_atoms_alpha = []
        ref_atoms_beta = []
        sample_atoms_alpha = []
        sample_atoms_beta = []
        first_ref = False
        second_ref = False
        first_sample = False
        second_sample = False

        for ref_chain in ref_model:
            if ref_chain.get_id() == chain_id_ref or chain_id_ref == 'ALL':
                if ref_chain.get_id() == chains_ref['ALPHA']:
                    for ref_res in ref_chain:
                            if not first_ref:
                                starting_count = ref_res.get_id()[1]
                                atoms_to_be_aligned = range(start_a + starting_count, end_a + starting_count - 1)
                                print('1: ')
                                print(atoms_to_be_aligned)
                                first_ref = True
                            if ref_res.get_id()[1] in atoms_to_be_aligned:
                                ref_atoms_alpha.append(ref_res['CA'])
                elif ref_chain.get_id() == chains_ref['BETA']:
                    for ref_res in ref_chain:
                            if not second_ref:
                                starting_count = ref_res.get_id()[1]
                                atoms_to_be_aligned = range(start_b + starting_count, end_b + starting_count - 1)
                                print('2: ')
                                print(atoms_to_be_aligned)
                                second_ref = True
                            if ref_res.get_id()[1] in atoms_to_be_aligned:
                                ref_atoms_beta.append(ref_res['CA'])

        for sample_chain in sample_model:
            if sample_chain.get_id() == chain_id_sample or chain_id_sample == 'ALL':
                if sample_chain.get_id() == chains_sample['ALPHA']:
                    for sample_res in sample_chain:
                        if not first_sample:
                            starting_count = sample_res.get_id()[1]
                            atoms_to_be_aligned = range(start_a + starting_count, end_a + starting_count - 1)
                            print('3: ')
                            print(atoms_to_be_aligned)
                            first_sample = True
                        if sample_res.get_id()[1] in atoms_to_be_aligned:
                            sample_atoms_alpha.append(sample_res['CA'])
                elif sample_chain.get_id() == chains_sample['BETA']:
                    for sample_res in sample_chain:
                        if not second_sample:
                            starting_count = sample_res.get_id()[1]
                            atoms_to_be_aligned = range(start_b + starting_count, end_b + starting_count - 1)
                            print('4: ')
                            print(atoms_to_be_aligned)
                            second_sample = True
                        if sample_res.get_id()[1] in atoms_to_be_aligned:
                            sample_atoms_beta.append(sample_res['CA'])

        ref_atoms = ref_atoms_alpha + ref_atoms_beta
        sample_atoms = sample_atoms_alpha + sample_atoms_beta

        print(len(ref_atoms))
        print(len(sample_atoms))

        super_imposer = PDB.Superimposer()
        super_imposer.set_atoms(ref_atoms, sample_atoms)
        super_imposer.apply(sample_atoms)

        # io = PDB.PDBIO()
        # io.set_structure(sample_structure)
        # io.save(self.get_pdb_id() + '_aligned.pdb')

        return super_imposer.rms

    def align_tcr_global(self, ref, sample):
        from Bio.pairwise2 import format_alignment
        matrix = matlist.blosum62
        self.set_file_name(ref)
        tcr_ref = self.get_tcr_amino_seq_V1('ALPHA') + '' + self.get_tcr_amino_seq_V1('BETA')
        self.set_file_name(sample)
        tcr_sample = self.get_tcr_amino_seq_V1('ALPHA') + '' + self.get_tcr_amino_seq_V1('BETA')
        # for a in pairwise2.align.globalmx(tcr_ref, tcr_sample, 2, -1):
        #     print(format_alignment(*a))
        score = pairwise2.align.globalmx(tcr_ref, tcr_sample, 2, -1, score_only=True)
        return score

    def get_gap(self, num):
        result = ''
        count = 0
        while count != num:
            result += '-'
            count += 1
        return result

    def cdr_align(self, ref, sample):
        self.set_file_name(ref)
        ref_seq = self.get_tcr_amino_seq_V1('ALPHA')
        ref_cdr1_alpha = self.get_cdr_seq(1, 'ALPHA')
        ref_cdr2_alpha = self.get_cdr_seq(2, 'ALPHA')
        ref_cdr1_beta = self.get_cdr_seq(1, 'BETA')
        ref_cdr2_beta = self.get_cdr_seq(2, 'BETA')
        ref_seq = ref_seq.replace(ref_cdr1_alpha, self.get_gap(len(ref_cdr1_alpha)))
        ref_seq = ref_seq.replace(ref_cdr2_alpha, self.get_gap(len(ref_cdr2_alpha)))
        print(ref_seq[0:138])
        print('a1 = ' + ref_cdr1_alpha + ' a2 = ' + ref_cdr2_alpha)
        self.set_file_name(sample)
        print(self.get_tcr_amino_seq_V1('ALPHA'))
        sample_cdr1_alpha = self.get_cdr_seq(1, 'ALPHA')
        sample_cdr2_alpha = self.get_cdr_seq(2, 'ALPHA')
        sample_cdr1_beta = self.get_cdr_seq(1, 'BETA')
        sample_cdr2_beta = self.get_cdr_seq(2, 'BETA')
        print('a1 = ' + sample_cdr1_alpha + ' a2 = ' + sample_cdr2_alpha)

    def pull_seq_pir(self, pdb_name):
        output = ''
        flag = False
        with open('ALL_TCR_PIR.pir') as f:
            for line in f:
                if line[4:].__contains__(pdb_name) and not flag:
                    flag = True
                elif line[0].__contains__('>') or line[0:2].__contains__('C;'):
                    flag = False
                elif flag is True and not line[0:9].__contains__('structure'):
                    output += line
        return output

    def prepare_seq(self, send):
        seq_a = self.get_tcr_amino_seq_V1('ALPHA')
        seq_b = self.get_tcr_amino_seq_V1('BETA')
        cdr1_a = self.get_cdr_seq(1, 'ALPHA')
        cdr2_a = self.get_cdr_seq(2, 'ALPHA')
        cdr1_b = self.get_cdr_seq(1, 'BETA')
        cdr2_b = self.get_cdr_seq(2, 'BETA')
        if send == 'a':
            return seq_a[0:115]
        elif send == 'b':
            return seq_b[0:115]
        elif send == '1a':
            return cdr1_a
        elif send == '2a':
            return cdr2_a
        elif send == '1b':
            return cdr1_b
        elif send == '2b':
            return cdr2_b
        elif send == 's1a':
            return seq_a[0: int(seq_a.find(cdr1_a))]
        elif send == 's2a':
            return seq_a[int(seq_a.find(cdr1_a)) + len(cdr1_a): int(seq_a.find(cdr2_a))]
        elif send == 's3a':
            count = int(seq_a.find(cdr2_a)) + len(cdr2_a)
            return seq_a[count:115]
        elif send == 's1b':
            return seq_b[0: int(seq_b.find(cdr1_b))]
        elif send == 's2b':
            return seq_b[int(seq_b.find(cdr1_b)) + len(cdr1_b): int(seq_b.find(cdr2_b))]
        elif send == 's3b':
            count = int(seq_b.find(cdr2_b)) + len(cdr2_b)
            return seq_b[count:115]

    def get_align_score(self, ref, sample):
        score = pairwise2.align.localms(ref, sample,  2, -1, -.5, -.1, score_only=True)
        return score

    # Returns a dictionary of highest scoring segments of TCR based on local alignment.
    # Allows for an exempt structure as a parameter.
    def prepare_align(self, exempt='****'):
        s1a = self.prepare_seq('s1a')
        s2a = self.prepare_seq('s2a')
        s3a = self.prepare_seq('s3a')
        s1b = self.prepare_seq('s1b')
        s2b = self.prepare_seq('s2b')
        s3b = self.prepare_seq('s3b')
        score_1 = []
        score_2 = []
        score_3 = []
        score_1b = []
        score_2b = []
        score_3b = []
        with open('PDB_NAMES.txt', 'r') as f:
            for name in f:
                if not name.__contains__(exempt):
                    sample_seq = self.pull_seq_pir(name).replace('\n', '')
                    test_score = self.get_align_score(sample_seq[:-2].replace('/', ''), s1a)
                    score_1.append([test_score, name[:-1]])
                    test_score = self.get_align_score(sample_seq[:-2].replace('/', ''), s2a)
                    score_2.append([test_score, name[:-1]])
                    test_score = self.get_align_score(sample_seq[:-2].replace('/', ''), s3a)
                    score_3.append([test_score, name[:-1]])
                    test_score = self.get_align_score(sample_seq[:-2].replace('/', ''), s1b)
                    score_1b.append([test_score, name[:-1]])
                    test_score = self.get_align_score(sample_seq[:-2].replace('/', ''), s2b)
                    score_2b.append([test_score, name[:-1]])
                    test_score = self.get_align_score(sample_seq[:-2].replace('/', ''), s3b)
                    score_3b.append([test_score, name[:-1]])
        sort_score = {'s1a': sorted(score_1[:-1])[-1][1], 's2a': sorted(score_2[:-1])[-1][1], 's3a': sorted(score_3[:-1])[-1][1], 's1b': sorted(score_1b[:-1])[-1][1], 's2b': sorted(score_2b[:-1])[-1][1], 's3b': sorted(score_3b[:-1])[-1][1]}
        return sort_score

    # Returns the specific sequences from each PDB that will be used for alignments.
    def prepare_info(self, exempt='****'):
        align = self.prepare_align(exempt)
        segments = {}
        for key in align:
            self.set_file_name('ALL_PDB/' + align[key] + '.pdb')
            if align[key] not in segments.keys():
                segments[align[key]] = [self.prepare_seq(key)]
            else:
                segments[align[key]].append(self.prepare_seq(key))
        return segments

    # Returns the alignment strings to run through alignment for each reference structure.
    def prepare_dash(self, exempt='****'):
        save_file = self.file_name
        info = self.prepare_info(exempt)
        print('4: ' + save_file)
        result = {}
        for key in info:
            self.set_file_name('ALL_PDB/' + key + '.pdb')
            alpha_chain = self.get_tcr_amino_seq_V1('ALPHA')
            new_alpha_chain = self.get_gap(len(alpha_chain))
            beta_chain = self.get_tcr_amino_seq_V1('BETA')
            new_beta_chain = self.get_gap(len(beta_chain))
            position = {}
            for seq in info[key]:
                if alpha_chain.__contains__(seq):
                    start = alpha_chain.find(seq)
                    end = len(seq) + start
                    if not position.__contains__('ALPHA'):
                        position['ALPHA'] = {seq: {'start': int(start), 'end': int(end)}}
                    else:
                        position['ALPHA'][seq] = {'start': int(start), 'end': int(end)}
                elif beta_chain.__contains__(seq):
                    start = beta_chain.find(seq)
                    end = len(seq) + start
                    if not position.__contains__('BETA'):
                        position['BETA'] = {seq: {'start': int(start), 'end': int(end)}}
                    else:
                        position['BETA'][seq] = {'start': int(start), 'end': int(end)}
            for chain in position:
                for seq in position[chain]:
                    if chain == 'ALPHA':
                        new_alpha_chain = new_alpha_chain[0:position[chain][seq]['start']] + seq \
                                          + new_alpha_chain[position[chain][seq]['end']:]
                    if chain == 'BETA':
                        new_beta_chain = new_beta_chain[0:position[chain][seq]['start']] + seq \
                                         + new_beta_chain[position[chain][seq]['end']:]
            result[key] = {'ALPHA': new_alpha_chain[0:115], 'BETA': new_beta_chain[0:115]}
        self.set_file_name(save_file)
        return result

    # Returns the number associated with the first atom in a chain
    def get_first_atom_on_chain(self, chain):
        with open(self.get_file_name(), 'r') as f:
            for line in f:
                if line[0:6] == 'ATOM  ':
                    if line[21] == chain:
                        return int(line[22:26])

    def create_ali_dash(self, dash_ali, chain='*****'):
        ori_pdb = self.get_file_name()
        print('3: ' + ori_pdb)
        output = ''
        for pdb in dash_ali:
            count = 1
            self.set_file_name('All_PDB/' + pdb + '.pdb')
            tcr_list = self.get_tcr_chains()
            if chain == 'ALPHA':
                first = self.get_first_atom_on_chain(tcr_list['ALPHA'])
                output += '>P1;' + self.get_pdb_id() + tcr_list['ALPHA'] + '\n'
                output += 'structureX:' + self.get_pdb_id() + ':' + str(first)
                output += ':' + tcr_list['ALPHA'] + ': ' + str(first + 115)
                output += ' :' + tcr_list['ALPHA'] + ': : :' + str(self.get_resolution()).rjust(5) + ":-1.00\n"
                for aa in dash_ali[pdb]['ALPHA']:
                    if count % 76 != 0:
                        output += aa
                        count += 1
                    else:
                        output += '\n' + aa
                        count += 1
                output += '*\n'
            # elif chain == 'BETA':
            #
            # else:
        self.set_file_name(ori_pdb)
        return output

    # Creates a new aliment file for the alpha chain based on the similar sequences found in prepare align
    # Uses an input PDB name for leave one out approach.
    def create_tcr_ali_alpha(self, exempt='****'):
        file_name = self.get_pdb_id() + '_alpha.ali'
        print(file_name)
        dash_ali = self.prepare_dash(exempt)
        print('1: ' + self.get_file_name())
        with open(file_name, 'w+') as f:
            print('2: ' + self.get_pdb_id())
            f.write(self.ali_alpha_text())
            f.write(self.create_ali_dash(dash_ali, 'ALPHA'))

    # Returns the alpha and beta chains with raw text for sequence to submit to modeller programs
    def raw_chains(self):
        tcr_alpha_chain = self.get_amino_acid_on_chain('A')
        tcr_beta_chain = self.get_amino_acid_on_chain('B')
        pdb_id = self.get_pdb_id()
        with open('STCRDab_Var_Raw', 'a+') as f:
            if tcr_alpha_chain != '':
                f.write('>' + pdb_id + 'A\n')
                for aa in tcr_alpha_chain:
                    f.write(aa)
                f.write('\n')
            if tcr_beta_chain != '':
                f.write('>' + pdb_id + 'B\n')
                for aa in tcr_beta_chain:
                    f.write(aa)
                f.write('\n')

    # Uses txt file returned from PISCES to compare resolution of structures and returns structures to remove.
    def redundant_resolution(self):
        remove = {}
        final_remove = []
        with open('STCRDab_9_1_abTCR_clean_LOG.txt', 'r') as log:
            for line in log:
                pdb1 = line[7:12].lower()
                pdb2 = line[13:18].lower()
                self.set_file_name('STCRDab_9_1_abTCR/' + pdb1[:-1] + '.pdb')
                resolution1 = self.get_resolution()
                self.set_file_name('STCRDab_9_1_abTCR/' + pdb2[:-1] + '.pdb')
                resolution2 = self.get_resolution()
                if resolution1 <= resolution2:
                    if pdb2[:-1] in remove:
                        if pdb1 not in remove[pdb2[:-1]]:
                            remove[pdb2[:-1]].append(pdb1)
                    else:
                        remove[pdb2[:-1]] = [pdb1]
                elif resolution1 == resolution2:
                    if pdb1[3] <= pdb2[3]:
                        if pdb2[:-1] in remove:
                            if pdb1 not in remove[pdb2[:-1]]:
                                remove[pdb2[:-1]].append(pdb1)
                        else:
                            remove[pdb2[:-1]] = [pdb1]
                    else:
                        if pdb1[:-1] in remove:
                            if pdb2 not in remove[pdb1[:-1]]:
                                remove[pdb1[:-1]].append(pdb2)
                        else:
                            remove[pdb1[:-1]] = [pdb2]
                else:
                    if pdb1[:-1] in remove:
                        if pdb2 not in remove[pdb1[:-1]]:
                            remove[pdb1[:-1]].append(pdb2)
                    else:
                        remove[pdb1[:-1]] = [pdb2]
        for key in remove:
            length = len(remove[key])
            for i in range(length):
                for k in range(length):
                    if k != i:
                        if remove[key][i][:-1] == remove[key][k][:-1]:
                            if key not in final_remove:
                                final_remove.append(key)
        # print(len(final_remove))
        # print(remove)
        return final_remove

    def create_test_folder(self, clean_folder, test_folder, cd_file):
        redundant_list = tool.cd_redundant_v2(cd_file)
        print(redundant_list)
        print(len(redundant_list))
        directory = '/Users/austinseamann/PycharmProjects/TCRDOCK/%s' % clean_folder
        copy_dir = '/Users/austinseamann/PycharmProjects/TCRDOCK/%s' % test_folder
        for filename in os.listdir(directory):
            if filename.endswith('.pdb'):
                if redundant_list.__contains__(filename[0:4]):
                    shutil.copyfile(directory + '/' + filename, copy_dir + '/' + filename)

    def create_blacklist_pisces(self, test_tcr):
        test_tcr_file = test_tcr + '.pdb'
        test_tcr = test_tcr.lower()
        file_path = 'STCRDab_ALL/%s' % test_tcr_file
        tool.set_file_name(file_path)
        tcr_list = tool.get_tcr_chains()
        alpha_list = []
        beta_list = []
        with open('STCRDab_ALL_Chains_LOG.txt', 'r+') as log:
            for line in log:
                pdb1 = line[7:12].lower()
                pdb2 = line[13:18].lower()
                if (test_tcr + tcr_list['ALPHA'].lower()) == pdb1:
                    alpha_list.append(pdb2[:-1])
                if (test_tcr + tcr_list['ALPHA'].lower()) == pdb2:
                    alpha_list.append(pdb1[:-1])
                if (test_tcr + tcr_list['BETA'].lower()) == pdb1:
                    beta_list.append(pdb2[:-1])
                if (test_tcr + tcr_list['BETA'].lower()) == pdb2:
                    beta_list.append(pdb1[:-1])
        both_list = set.intersection(set(alpha_list), set(beta_list))
        return list(both_list)

    def trim_tcr(self):
        matrix = matlist.blosum62
        inds = {'ALPHA': [], 'BETA': []}
        tcr_list = self.get_tcr_chains()
        chains = {'ALPHA': self.get_amino_acid_on_chain(tcr_list['ALPHA']),
                  'BETA': self.get_amino_acid_on_chain(tcr_list['BETA'])}  # Pulls aa seq. for alpha/beta chains
        cuts = {'ALPHA': ['QFGAG', 'VTP'], 'BETA': ['FGPGT', 'VTE']}  # Cut alignment (first is normal) trying to implement cutoff at start of chain.
        control_chains = {'ALPHA': [
            'KEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFTYRQYSGKSPELIMSIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSDSATYLCAVTTDSTGKLQFGAGTQVVVTPDIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVATSNKSDFACANAFNNSIIPEDTFFPSPESS',
            'QKVTQTQTSISVMEKTTVTMDCVYETQDSSYFLFTYKQTASGEIVFLIRQDSYKKENATVGHYSLNFQKPKSSIGLIITATQIEDSAVYFCAMRGDYGGSGNKLIFGTGTLLSVKP'],
            'BETA': [
            'NAGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSTYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVTEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSTTVNGKEVHSGVSTDPQPLKEQPALNDSRYALSSRLRVSATFTQNPRNHFRCQVQFYGLSENDETTQDRAKPVTQIVSAEATGRAD',
            'VTLLEQNPRTRLVPRGQAVNLRCILKNSQYPTMSTYQQDLQKQLQTLFTLRSPGDKEVKSLPGADYLATRVTDTELRLQVANMSQGRTLYCTCSADRVGNTLYFGEGSRLIV']}
        for chain in chains:
            for cut in cuts[chain]:
                lc = len(cut)
                window, count = self.big_frame(control_chains[chain][0], chains[chain], cut)
                tmp = []
                for s in range(len(window) - lc):
                    score = pairwise2.align.globaldx(cut, window[s:s+lc], matrix, score_only=True)
                    seq2 = difflib.SequenceMatcher(None, cut, window[s:s+lc])
                    perc = seq2.ratio() * 100
                    tmp.append([score, perc, count, window[s:s+lc]])  # resulting lists
                inds[chain].append(max(tmp))  # Where it's freezing up right now.
        return inds

    def big_frame(self, control, seq, c):
        matrix = matlist.blosum62
        ind0 = control.find(c)
        align = pairwise2.align.globaldx(control, seq, matrix, one_alignment_only=True)
        control_align = align[0][0]
        sample_align = align[0][1]
        count = 0
        count2 = 0
        for i in control_align:
            if count < ind0:
                if i != '-':
                    count += 1
                    count2 += 1
                else:
                    count2 += 1
            else:
                break
        window = sample_align[count2 - 25:count2 + 25]
        window_seq = window.replace('-', '')
        return window_seq, count

    # Takes input of two .clstr file from CD-Hit result and creates save_list based on higher resolution
    def cd_redundant_v1(self, cluster_a, cluster_b):
        alpha = {}
        beta = {}
        alpha_best = {}
        beta_best = {}
        alpha_save = set()
        beta_save = set()
        save_list = []
        with open(cluster_a, 'r') as file_a:
            for line in file_a:
                if line[0] == '>':
                    cluster = line[1:-1]
                    alpha[cluster] = []
                else:
                    index = line.index('>') + 1
                    alpha[cluster].append(line[index:index+4])
        with open(cluster_b, 'r') as file_b:
            for line in file_b:
                if line[0] == '>':
                    cluster = line[1:-1]
                    beta[cluster] = []
                else:
                    index = line.index('>') + 1
                    beta[cluster].append(line[index:index+4])
        for key in alpha:
            res_dict = {}
            if len(alpha[key]) == 1:
                alpha_save.add(alpha[key][0])
            else:
                for pdb in alpha[key]:
                    self.set_file_name('STCRDab_ALL/' + pdb + '.pdb')
                    res = self.get_resolution()
                    res_dict[pdb] = res
                sort_res_dict = sorted(res_dict.items(), key=lambda x: x[1])
                best = sort_res_dict[0][0]
                alpha_best[best] = []
                for each in sort_res_dict:
                    alpha_best[best].append(each[0])
        for key in beta:
            res_dict = {}
            if len(beta[key]) == 1:
                beta_save.add(beta[key][0])
            else:
                for pdb in beta[key]:
                    self.set_file_name('STCRDab_ALL/' + pdb + '.pdb')
                    res = self.get_resolution()
                    res_dict[pdb] = res
                sort_res_dict = sorted(res_dict.items(), key=lambda x: x[1])
                best = sort_res_dict[0][0]
                beta_best[best] = []
                for each in sort_res_dict:
                    beta_best[best].append(each[0])
        for key in alpha_best:
            if key in beta_best.keys():
                save_list.append(key)
            # elif key in beta_best.items():
            #     save_list.append(key)
        union = alpha_save.union(beta_save)
        save_list.append(union)
        print(alpha)
        print(beta)
        print(union)
        print(len(union))
        print(alpha_best)
        print(beta_best)
        print(save_list)
        print(len(save_list))

    # Returns a list of TCRs with higher resolution based on inputted clstr file.
    def cd_redundant_v2(self, cluster):
        tcr = {}
        tcr_best = {}
        save_list = []
        with open(cluster, 'r') as file:
            for line in file:
                if line[0] == '>':
                    cluster = line[1:-1]
                    tcr[cluster] = []
                else:
                    index = line.index('>') + 1
                    tcr[cluster].append(line[index:index+4])
        for key in tcr:
            res_dict = {}
            if len(tcr[key]) == 1:
                save_list.append(tcr[key][0])
            else:
                for pdb in tcr[key]:
                    self.set_file_name('ALL_PDB/' + pdb + '.pdb')
                    res = self.get_resolution()
                    res_dict[pdb] = res
                sort_res_dict = sorted(res_dict.items(), key=lambda x: x[1])
                best = sort_res_dict[0][0]
                tcr_best[best] = []
                for each in sort_res_dict:
                    tcr_best[best].append(each[0])
        for key in tcr_best:
            save_list.append(key)
        return save_list

    def cd_blacklist(self, cluster_a, cluster_b):
        alpha = {}
        beta = {}
        black_list = {}
        pdb = self.get_pdb_id()
        with open(cluster_a, 'r') as file_a:
            for line in file_a:
                if line[0] == '>':
                    cluster = line[1:-1]
                    alpha[cluster] = []
                else:
                    index = line.index('>') + 1
                    alpha[cluster].append(line[index:index + 4])
                    if line[index:index + 4] == pdb:
                        alpha_cluster = cluster
        with open(cluster_b, 'r') as file_b:
            for line in file_b:
                if line[0] == '>':
                    cluster = line[1:-1]
                    beta[cluster] = []
                else:
                    index = line.index('>') + 1
                    beta[cluster].append(line[index:index + 4])
                    if line[index:index + 4] == pdb:
                        beta_cluster = cluster
        black_list['ALPHA'] = alpha[alpha_cluster]
        black_list['BETA'] = beta[beta_cluster]
        return black_list

    # Setup right now to produce Fasta files for the entire chain and not the trimmed structures.
    def rep_builder_entry(self, cluster_a, cluster_b, percent):
        directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Testing_full"
        if percent == 90:
            alpha_file = 'alpha_90_full.fasta'
            beta_file = 'beta_90_full.fasta'
        elif percent == 80:
            alpha_file = 'alpha_80_full.fasta'
            beta_file = 'beta_80_full.fasta'
        with open('Modeled/RepBuilder/' + alpha_file, 'w+') as file_a:
            with open('Modeled/RepBuilder/' + beta_file, 'w+') as file_b:
                for filename in os.listdir(directory):
                    if filename.endswith(".pdb"):
                        self.set_file_name('STCRDab_abTCR_clean/' + filename)
                        tcr_alpha_chain = self.get_amino_acid_on_chain('A')
                        tcr_beta_chain = self.get_amino_acid_on_chain('B')
                        pdb_id = self.get_pdb_id()
                        black_list = self.cd_blacklist(cluster_a, cluster_b)
                        set_a = set(black_list['ALPHA'])
                        set_b = set(black_list['BETA'])
                        union = set_a.union(set_b)
                        union_string = ''
                        for each in union:
                            union_string += each + ','
                        count_1 = 1
                        count_2 = 1
                        if tcr_alpha_chain != '':
                            file_a.write('>' + pdb_id + ' ((' + union_string[:-1] + '))\n')
                            for aa in tcr_alpha_chain:
                                if count_1 % 81 != 0:
                                    file_a.write(aa)
                                    count_1 += 1
                                else:
                                    file_a.write('\n' + aa)
                                    count_1 = 2
                            file_a.write('\n')
                        if tcr_beta_chain != '':
                            file_b.write('>' + pdb_id + ' ((' + union_string[:-1] + '))\n')
                            for aa in tcr_beta_chain:
                                if count_2 % 81 != 0:
                                    file_b.write(aa)
                                    count_2 += 1
                                else:
                                    file_b.write('\n' + aa)
                                    count_2 = 2
                            file_b.write('\n')

    def other_entry(self, cluster_a, cluster_b, percent):
        directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Testing_full"
        sort_dir = sorted(os.listdir(directory))
        if percent == 90:
            txt_file = 'Modeled/all_90.txt'
        elif percent == 80:
            txt_file = 'Modeled/all_80.txt'
        with open(txt_file, 'w+') as file:
            for filename in sort_dir:
                if filename.endswith(".pdb"):
                    self.set_file_name('STCRDab_abTCR_clean/' + filename)
                    tcr_alpha_chain = self.get_amino_acid_on_chain('A')
                    tcr_beta_chain = self.get_amino_acid_on_chain('B')
                    pdb_id = self.get_pdb_id()
                    black_list = self.cd_blacklist(cluster_a, cluster_b)
                    set_a = set(black_list['ALPHA'])
                    set_b = set(black_list['BETA'])
                    union = set_a.union(set_b)
                    union_string = ''
                    for each in union:
                        union_string += each + ','
                    file.write(pdb_id + '\n' + 'ALPHA:\n' + tcr_alpha_chain + '\nBETA:\n' + tcr_beta_chain + '\n')
                    file.write('Black List:\n' + union_string[:-1] + '\n\n')

    # LYRA/TCRmodel = A/B RepBuilder = L/H
    def check_seq(self, modeled_dir):
        test_seq = {}
        test_seq['ALPHA'] = self.get_amino_acid_on_chain('A')
        test_seq['BETA'] = self.get_amino_acid_on_chain('B')
        pdb_id = self.get_pdb_id()
        file_name = pdb_id + '.pdb'
        if file_name in os.listdir('/Users/austinseamann/PycharmProjects/TCRDOCK/' + modeled_dir):
            print(pdb_id)
            self.set_file_name(modeled_dir + '/' + file_name)
            modeled_seq = {}
            modeled_seq['ALPHA'] = self.get_amino_acid_on_chain('A')
            modeled_seq['BETA'] = self.get_amino_acid_on_chain('B')
            print('ALPHA')
            print('CRY:\n' + test_seq['ALPHA'])
            for a in pairwise2.align.globalms(test_seq['ALPHA'], modeled_seq['ALPHA'], 2, -1, -2, -.5, penalize_end_gaps=(False, False), one_alignment_only=True):
                print(pairwise2.format_alignment(*a))
            print('BETA')
            print('CRY:\n' + test_seq['BETA'])
            for a in pairwise2.align.globalms(test_seq['BETA'], modeled_seq['BETA'], 2, -1, -2, -.5, penalize_end_gaps=(False, False), one_alignment_only=True):
                print(pairwise2.format_alignment(*a))

    # Initialize the program using the modeled PDBs, will produce a directory for newly trimmed crystal structures
    # using a semi-global alignment to remove the start and ends of structure that wasn't modeled by muting the aa.
    error_list = []  # Gloabl variable to produce an error list for semi_align_pdb
    def semi_align_pdb(self, pdb):
        modeled_alpha = self.get_amino_acid_on_chain('A')
        modeled_beta = self.get_amino_acid_on_chain('B')
        pdb_id = pdb
        print(pdb_id)
        self.set_file_name('Testing_align/' + pdb_id + '_tcr.pdb')
        cry_alpha = self.get_amino_acid_on_chain('A')
        cry_beta = self.get_amino_acid_on_chain('B')
        alpha_align = pairwise2.align.globalms(cry_alpha, modeled_alpha, 2, -1, -2, -.5, penalize_end_gaps=(False, False), one_alignment_only=True)
        beta_align = pairwise2.align.globalms(cry_beta, modeled_beta, 2, -1, -2, -.5, penalize_end_gaps=(False, False), one_alignment_only=True)
        try:
            left_a, right_a = self.gap_finder(alpha_align[0][1])
            left_b, right_b = self.gap_finder(beta_align[0][1])
            self.mute_aa(0, left_a, 'A')
            print('ALPHA-Gap: ' + alpha_align[0][1])
            print('left_a: ' + str(left_a))
            self.mute_aa(right_a, 1000, 'A')
            print('right_a: ' + str(right_a))
            print('BETA-Gap: ' + beta_align[0][1])
            self.mute_aa(0, left_b, 'B')
            print('left_b:' + str(left_b))
            self.mute_aa(right_b, 1000, 'B')
            print('right_b: ' + str(right_b))
            if right_a == 0 or right_b == 0:
                error_list.append(pdb_id)
        except TypeError:
            print('TypeError')
            error_list.append(pdb_id)

    def unmute_all(self, dir):
        for filename in os.listdir(dir):
            if filename.endswith(".pdb"):
                with open(dir + filename, 'r') as r:
                    data = r.readlines()
                with open(dir + filename, 'w+') as w:
                    for line in data:
                        if line[0:6] == 'DEATOM':
                            result = line.replace('DEATOM', 'ATOM  ', 1)
                            w.write(result)
                        else:
                            w.write(line)

    # Input of an aligned sequence in, tells the start and end positions of the aligned portion and returns if gap in
    # the middle of the alignment.
    def gap_finder(self, seq_in):
        flag_1 = False
        flag_2 = False
        start_count = 0
        end_count = 0
        for letter in seq_in:
            if letter == '-' and not flag_1 and not flag_2:  # Regular start
                start_count += 1
                end_count += 1
                flag_1 == True
            elif letter == '-' and flag_1 and not flag_2:  # After regular start
                start_count += 1
                end_count += 1
            elif letter != '-' and flag_1 and not flag_2:  # Once we hit our first aa
                end_count += 1
                flag_2 = True
            elif letter != '-' and flag_1 and flag_2:  # Continuing aa count
                end_count += 1
            elif letter == '-' and flag_1 and flag_2:  # End of seq
                #  Removed end_count += 1
                dif = end_count - start_count
                if dif <= 100:
                    print('Failed: ' + self.get_pdb_id())
                    return 0, 0
                else:
                    return start_count, end_count
            elif letter != '-' and not flag_1 and not flag_2:  # Start of seq is 0
                flag_1 = True
                flag_2 = True
                end_count += 1
        return start_count, end_count

    # Mutes atoms of amino acid positions to trim tcr based on left, right, and chain_id
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

    # Unmutes atoms of amino acid positions to untrim tcr based on left, right, and chain_id
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

    # Switches around the beta and alpha chain in PDB file for input into RMSD program.
    def clean_repbuilder(self):
        beta_lines = ''
        alpha_lines = ''
        header = ''
        result = ''
        atom_count = 0
        print(self.get_file_name())
        with open(self.get_file_name(), 'r') as read:
            for line in read:
                if line[0:6] == 'ATOM  ' or line[0:6] == 'TER   ':
                    if line[21] == 'B':
                        beta_lines += line
                    if line[21] == 'A':
                        alpha_lines += line
                else:
                    header += line
        alpha = iter(alpha_lines.splitlines())
        for line in alpha:
            atom_count += 1
            num = line[6:11]
            result += line.replace(num, str(atom_count).rjust(5), 1)
            result += '\n'
        beta = iter(beta_lines.splitlines())
        for line in beta:
            atom_count += 1
            num = line[6:11]
            result += line.replace(num, str(atom_count).rjust(5), 1)
            result += '\n'
        with open(self.get_file_name(), 'w+') as w:
            w.write(header)
            w.write(result)

    # Returns longest common subsequence of three strings
    def lcs3(self, a, b, c):
        m = len(a)
        l = len(b)
        n = len(c)
        subs = [[[0 for k in range(n + 1)] for j in range(l + 1)] for i in range(m + 1)]

        for i, x in enumerate(a):
            for j, y in enumerate(b):
                for k, z in enumerate(c):
                    if x == y and y == z:
                        subs[i + 1][j + 1][k + 1] = subs[i][j][k] + 1
                    else:
                        subs[i + 1][j + 1][k + 1] = max(subs[i + 1][j + 1][k],
                                                        subs[i][j + 1][k + 1],
                                                        subs[i + 1][j][k + 1])
        # return subs[-1][-1][-1] #if you only need the length of the lcs
        lcs = ""
        while m > 0 and l > 0 and n > 0:
            step = subs[m][l][n]
            if step == subs[m - 1][l][n]:
                m -= 1
            elif step == subs[m][l - 1][n]:
                l -= 1
            elif step == subs[m][l][n - 1]:
                n -= 1
            else:
                lcs += str(a[m - 1])
                m -= 1
                l -= 1
                n -= 1

        return lcs[::-1]

    # Determines longest common sequence between LYRA/RepBuilder/TCRModel models.
    # Assumes that models are in "Modeled/*/--" and reference structures will be saved to "Testing_align"
    # -- = %identity and * will do all modeling pipelines
    def lcs_pdbs(self, percentage):
        lyra_dir = "Modeled/LYRA/" + str(percentage) + '/'
        repbuilder_dir = "Modeled/RepBuilder/" + str(percentage) + '/'
        tcrmodel_dir = "Modeled/TCRModel/" + str(percentage) + '/'
        lyra_pdb = os.listdir(lyra_dir)
        lyra_pdb_list = set([fname for fname in lyra_pdb if fname.endswith('.pdb')])
        repbuilder_pdb = os.listdir(repbuilder_dir)
        repbuilder_pdb_list = set([fname for fname in repbuilder_pdb if fname.endswith('.pdb')])
        tcrmodel_pdb = os.listdir(tcrmodel_dir)
        tcrmodel_pdb_list = set([fname for fname in tcrmodel_pdb if fname.endswith('.pdb')])
        intersection_list = set.intersection(lyra_pdb_list, repbuilder_pdb_list, tcrmodel_pdb_list)  # Intersection of all three pipeline outputs
        print(intersection_list)
        lyra_dic = {}
        repbuilder_dic = {}
        tcrmodel_dic = {}
        lcs_dic = {}
        for pdb in sorted(intersection_list):
            print(pdb)
            self.set_file_name(lyra_dir + pdb)
            lyra_dic[pdb] = {'ALPHA': self.get_amino_acid_on_chain('A'), 'BETA': self.get_amino_acid_on_chain('B')}
            self.set_file_name(repbuilder_dir + pdb)
            repbuilder_dic[pdb] = {'ALPHA': self.get_amino_acid_on_chain('A'),
                                   'BETA': self.get_amino_acid_on_chain('B')}
            self.set_file_name(tcrmodel_dir + pdb)
            tcrmodel_dic[pdb] = {'ALPHA': self.get_amino_acid_on_chain('A'), 'BETA': self.get_amino_acid_on_chain('B')}
            lcs_alpha = self.lcs3(lyra_dic[pdb]['ALPHA'], repbuilder_dic[pdb]['ALPHA'], tcrmodel_dic[pdb]['ALPHA'])
            lcs_beta = self.lcs3(lyra_dic[pdb]['BETA'], repbuilder_dic[pdb]['BETA'], tcrmodel_dic[pdb]['BETA'])
            lcs_dic[pdb] = {'ALPHA': lcs_alpha, 'BETA': lcs_beta}
        with open('lcs_result_80.txt', 'w') as file:
            for pdb in lcs_dic:
                file.write(pdb + '\n')
                file.write('ALPHA ' + lcs_dic[pdb]['ALPHA'] + '\n')
                file.write('BETA ' + lcs_dic[pdb]['BETA'] + '\n')

    # Does an alignment on both the Crystal and Modeled structures based on the lcs of the three pipeline outputs
    # Expects input file to be 'pdb name\n' + 'ALPHA Seq\n' + 'BETA Seq\n'
    def lcs_semi_align(self, lcs_file, percentage, pipeline):
        lcs_dic = {}
        modeled_dir = 'Modeled_align/' + pipeline + "/" + str(percentage) + "/"
        with open(lcs_file, 'r') as file:
            for line in file:
                if line.endswith('.pdb\n'):
                    lcs_dic[line[:-1]] = {}
                    pdb_name = line[:-1]
                elif line[:5] == 'ALPHA':
                    lcs_dic[pdb_name]['ALPHA'] = line[6:-1]
                elif line[:4] == 'BETA':
                    lcs_dic[pdb_name]['BETA'] = line[5:-1]
        print(lcs_dic)
        for pdb in lcs_dic:
            lcs_alpha = lcs_dic[pdb]['ALPHA']
            lcs_beta = lcs_dic[pdb]['BETA']
            self.set_file_name(modeled_dir + pdb)
            modeled_alpha = self.get_amino_acid_on_chain('A')
            modeled_beta = self.get_amino_acid_on_chain('B')
            print(pdb + ' modeled')
            self.align_mute(modeled_alpha, lcs_alpha, modeled_beta, lcs_beta)
            self.set_file_name('Testing_align/' + pdb[:-4] + '_tcr.pdb')
            cry_alpha = self.get_amino_acid_on_chain('A')
            cry_beta = self.get_amino_acid_on_chain('B')
            print(pdb + ' cry')
            self.align_mute(cry_alpha, lcs_alpha, cry_beta, lcs_beta)

    def align_mute(self, muting_alpha, seq_alpha, muting_beta, seq_beta):
        print('Muting ALPHA: ' + muting_alpha)
        print('Seq ALPHA: ' + seq_alpha)
        print('Muting BETA: ' + muting_beta)
        print('Seq BETA: ' + seq_beta)
        alpha_align = pairwise2.align.globalms(muting_alpha, seq_alpha, 2, -1, -2, -.5, penalize_end_gaps=(False, False), one_alignment_only=True)
        beta_align = pairwise2.align.globalms(muting_beta, seq_beta, 2, -1, -2, -.5, penalize_end_gaps=(False, False), one_alignment_only=True)
        try:
            left_a, right_a = self.gap_finder(alpha_align[0][1])
            left_b, right_b = self.gap_finder(beta_align[0][1])
            self.mute_aa(0, left_a, 'A')
            print('ALPHA-Gap: ' + alpha_align[0][1])
            print('left_a: ' + str(left_a))
            if right_a != len(muting_alpha):
                self.mute_aa(right_a, 1000, 'A')
            print('right_a: ' + str(right_a))
            print('BETA-Gap: ' + beta_align[0][1])
            self.mute_aa(0, left_b, 'B')
            print('left_b:' + str(left_b))
            if right_b != len(muting_beta):
                self.mute_aa(right_b, 1000, 'B')
            print('right_b: ' + str(right_b))
            if right_a == 0 or right_b == 0:
                error_list.append(self.get_file_name())
        except TypeError:
            print('TypeError')
            error_list.append(self.get_file_name())
        except:
            print('OtherError')
            error_list.append(self.get_file_name())

    def remove_dot(self, input_string):
        output = ""
        for letter in input_string:
            if letter != ".":
                output += letter
        return output

    # Returns a dict. of all CDRs of TCRs available in IMGT database.
    def make_cdr_ref(self):
        pdb_cdr_dic = {}
        with open('IMGT_CDR.txt', 'r') as file:
            for line in file:
                listed = line.split()
                # print(listed)
                try:
                    value = int(listed[0])
                    if not pdb_cdr_dic.__contains__(listed[1][:-2]):
                        pdb_cdr_dic[listed[1][:-2]] = {'ALPHA': {'CDR1': None, 'CDR2': None, 'CDR3': None}, 'BETA': {'CDR1': None, 'CDR2': None, 'CDR3': None}}
                    for each in range(len(listed)):
                        if listed[each] == "V-ALPHA":
                            pdb_cdr_dic[listed[1][:-2]]['ALPHA']['CDR1'] = self.remove_dot(listed[each + 5])
                            pdb_cdr_dic[listed[1][:-2]]['ALPHA']['CDR2'] = self.remove_dot(listed[each + 8])
                            pdb_cdr_dic[listed[1][:-2]]['ALPHA']['CDR3'] = self.remove_dot(listed[each + 13])
                            break
                        elif listed[each] == "V-BETA":
                            pdb_cdr_dic[listed[1][:-2]]['BETA']['CDR1'] = self.remove_dot(listed[each + 5])
                            pdb_cdr_dic[listed[1][:-2]]['BETA']['CDR2'] = self.remove_dot(listed[each + 8])
                            pdb_cdr_dic[listed[1][:-2]]['BETA']['CDR3'] = self.remove_dot(listed[each + 13])
                            break
                except:
                    pass
        return pdb_cdr_dic

    # Returns the positions of the amino acids in the alignment chain.
    def gap_finder_cdr(self, seq_in):
        flag_1 = False
        flag_2 = False
        start_count = 0
        end_count = 0
        for letter in seq_in:
            if letter == '-' and not flag_1 and not flag_2:  # Regular start
                start_count += 1
                end_count += 1
                flag_1 == True
            elif letter == '-' and flag_1 and not flag_2:  # After regular start
                start_count += 1
                end_count += 1
            elif letter != '-' and flag_1 and not flag_2:  # Once we hit our first aa
                end_count += 1
                flag_2 = True
            elif letter != '-' and flag_1 and flag_2:  # Continuing aa count
                end_count += 1
            elif letter == '-' and flag_1 and flag_2:  # End of seq
                #  Removed end_count += 1
                dif = end_count - start_count
                if dif <= 2:
                    print('Failed: ' + self.get_pdb_id())
                    return 0, 0
            elif letter != '-' and not flag_1 and not flag_2:  # Start of seq is 0
                flag_1 = True
                flag_2 = True
                end_count += 1
        return start_count, end_count

    # Local alignment to find CDR sequences
    def local_align(self, large, small):
        matrix = matlist.blosum62
        align = pairwise2.align.localms(large, small, 1, -2, -2, 0, one_alignment_only=True)
        return self.gap_finder_cdr(align[0][1])

    # Finds position of CDR per cdr_type and chain_type based on relabeled PDB: Returns start, end
    def find_cdr_position(self, cdr_type, chain_type, cdr_list, pdb_id):
        if chain_type == 'ALPHA':
            seq = cdr_list[pdb_id]['ALPHA'][cdr_type]
            start, end = self.local_align(self.get_amino_acid_on_chain('A'), seq)
        elif chain_type == 'BETA':
            seq = cdr_list[pdb_id]['BETA'][cdr_type]
            start, end = self.local_align(self.get_amino_acid_on_chain('B'), seq)
        return start, end

    def create_zone_entry(self):
        file_name = self.get_file_name()[-12:]
        line_list = []
        LYRA_90 = "Modeled_align/LYRA/90/"
        LYRA_80 = "Modeled_align/LYRA/80/"
        Rep_90 = "Modeled_align/RepBuilder/90/"
        Rep_80 = "Modeled_align/RepBuilder/80/"
        TCR_90 = "Modeled_align/TCRModel/90/"
        TCR_80 = "Modeled_align/TCRModel/80/"
        dir_list = [LYRA_90, LYRA_80, Rep_90, Rep_80, TCR_90, TCR_80]
        pdb_name = file_name[:4]
        ccdr1_a = self.find_cdr_position('CDR1', 'ALPHA', cdr_list, pdb_name)
        ccdr2_a = self.find_cdr_position('CDR2', 'ALPHA', cdr_list, pdb_name)
        ccdr3_a = self.find_cdr_position('CDR3', 'ALPHA', cdr_list, pdb_name)
        ccdr1_b = self.find_cdr_position('CDR1', 'BETA', cdr_list, pdb_name)
        ccdr2_b = self.find_cdr_position('CDR2', 'BETA', cdr_list, pdb_name)
        ccdr3_b = self.find_cdr_position('CDR3', 'BETA', cdr_list, pdb_name)
        for dir in dir_list:
            try:
                self.set_file_name(dir + pdb_name + '.pdb')
                cdr1_a = self.find_cdr_position('CDR1', 'ALPHA', cdr_list, pdb_name)
                cdr2_a = self.find_cdr_position('CDR2', 'ALPHA', cdr_list, pdb_name)
                cdr3_a = self.find_cdr_position('CDR3', 'ALPHA', cdr_list, pdb_name)
                cdr1_b = self.find_cdr_position('CDR1', 'BETA', cdr_list, pdb_name)
                cdr2_b = self.find_cdr_position('CDR2', 'BETA', cdr_list, pdb_name)
                cdr3_b = self.find_cdr_position('CDR3', 'BETA', cdr_list, pdb_name)
                zone1 = 'A' + str(ccdr1_a[0]) + '-A' + str(ccdr1_a[1]) + ':A' + str(cdr1_a[0]) + '-A' + str(cdr1_a[1])
                zone2 = 'A' + str(ccdr2_a[0]) + '-A' + str(ccdr2_a[1]) + ':A' + str(cdr2_a[0]) + '-A' + str(cdr2_a[1])
                zone3 = 'A' + str(ccdr3_a[0]) + '-A' + str(ccdr3_a[1]) + ':A' + str(cdr3_a[0]) + '-A' + str(cdr3_a[1])
                zone4 = 'B' + str(ccdr1_b[0]) + '-B' + str(ccdr1_b[1]) + ':B' + str(cdr1_b[0]) + '-B' + str(cdr1_b[1])
                zone5 = 'B' + str(ccdr2_b[0]) + '-B' + str(ccdr2_b[1]) + ':B' + str(cdr2_b[0]) + '-B' + str(cdr2_b[1])
                zone6 = 'B' + str(ccdr3_b[0]) + '-B' + str(ccdr3_b[1]) + ':B' + str(cdr3_b[0]) + '-B' + str(cdr3_b[1])
                line_list.extend([zone1, zone2, zone3, zone4, zone5, zone6])
            except OSError:
                line_list.extend(['bad', 'bad', 'bad', 'bad', 'bad', 'bad'])
        with open('ProFit_CDR.txt', 'a+') as file:
            file.write(pdb_name + ' ')
            for zone in line_list:
                file.write(zone + ' ')
            file.write('\n')


tool = PdbTools()

# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Testing_full"
# cdr_list = tool.make_cdr_ref()
# for each in sorted(os.listdir(directory)):
#     if each.endswith(".pdb"):
#         tool.set_file_name(directory + "/" + each)
#         tool.create_zone_entry()
        # print(each)
        # print('CDR1_A')
        # print(tool.find_cdr_position('CDR1', 'ALPHA', cdr_list))
        # print('CDR2_A')
        # print(tool.find_cdr_position('CDR2', 'ALPHA', cdr_list))
        # print('CDR3_A')
        # print(tool.find_cdr_position('CDR3', 'ALPHA', cdr_list))
        # print('CDR1_B')
        # print(tool.find_cdr_position('CDR1', 'BETA', cdr_list))
        # print('CDR2_B')
        # print(tool.find_cdr_position('CDR2', 'BETA', cdr_list))
        # print('CDR3_B')
        # print(tool.find_cdr_position('CDR3', 'BETA', cdr_list))

# Step 1
# Going from STCRDab alpha/beta TCRs to cleaned structures
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/STCRDab_10_10_abTCR"
# for filename in os.listdir(directory):
#     if filename.endswith(".pdb"):
#         tool.set_file_name('STCRDab_10_10_abTCR/' + filename)
#         print(filename)
#         tool.clean_tcr_count('STCRDab_abTCR_clean/')

# Step 2 from clean structures to trimmed structures to submit to CD hit
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/STCRDab_abTCR_clean"
# for filename in os.listdir(directory):
#     if filename.endswith(".pdb"):
#         tool.set_file_name('STCRDab_abTCR_clean/' + filename)
#         print(filename)
#         tool.clean_tcr_count_trim('STCRDab_Trim/')

# Step 3
# Use for creating fasta files for submission to CD-Hit
# Activate conda environment "conda activate TCRDOCK"
# Run "cd-hit -i 99TCR.fasta -o 99TCR_CD -c 0.99 -n 5 -M 8000 -d 0" after
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/STCRDab_Trim"
# for filename in os.listdir(directory):
#     if filename.endswith(".pdb"):
#         tool.set_file_name('STCRDab_Trim/' + filename)
#         print(filename)
#         tool.fasta_TCR('CD-Hit/99/99TCR.fasta')

# Step 4
# Takes cd-hit clstr file and runs through cd_redundant_v2 method to create a folder with testing TCRs
# Change entry folder to STCRDab_Var for just variable region and STCRDab_9_1_abTCR_clean for full.
# tool.create_test_folder('STCRDab_abTCR_clean', 'Testing_full', 'CD-Hit/99/99TCR_CD.clstr')

# Step 5
# Use to write fasta files for submitting to CD-Hit for separate Alpha/Beta files
# Run "cd-hit -i TCR_Alpha.fasta -o TCR_Alpha_90 -c 0.90 -n 5 -M 8000 -d 0" after
# And "cd-hit -i TCR_Beta.fasta -o TCR_Beta_90 -c 0.90 -n 5 -M 8000 -d 0"
# 80: "cd-hit -i TCR_Alpha.fasta -o TCR_Alpha_80 -c 0.80 -n 5 -M 8000 -d 0"
# 80: "cd-hit -i TCR_Beta.fasta -o TCR_Beta_80 -c 0.80 -n 5 -M 8000 -d 0"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/STCRDab_Trim"
# for filename in os.listdir(directory):
#     if filename.endswith(".pdb"):
#         tool.set_file_name('STCRDab_Trim/' + filename)
#         print(filename)
#         tool.fasta_TCR('CD-Hit/TCR_Alpha.fasta', 'ALPHA')
#         tool.fasta_TCR('CD-Hit/TCR_Beta.fasta', 'BETA')

# Step 6
# Creating blacklist input for RB
# tool.rep_builder_entry('CD-Hit/TCR_Alpha_90.clstr', 'CD-Hit/TCR_Beta_90.clstr', 90)
# tool.other_entry('CD-Hit/TCR_Alpha_90.clstr', 'CD-Hit/TCR_Beta_90.clstr', 90)
# tool.rep_builder_entry('CD-Hit/TCR_Alpha_80.clstr', 'CD-Hit/TCR_Beta_80.clstr', 80)
# tool.other_entry('CD-Hit/TCR_Alpha_80.clstr', 'CD-Hit/TCR_Beta_80.clstr', 80)

# Step 6.5 - Run through modeling programs.

# Step 7
# Just swapping the Beta and Alpha chain order in RepBuilder results
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled/RepBuilder/90"
# for filename in os.listdir(directory):
#     if filename.endswith('.pdb'):
#         tool.set_file_name('Modeled/RepBuilder/90/' + filename)
#         tool.clean_repbuilder()
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled/RepBuilder/80"
# for filename in os.listdir(directory):
#     if filename.endswith('.pdb'):
#         tool.set_file_name('Modeled/RepBuilder/80/' + filename)
#         tool.clean_repbuilder()

# Copy over to Modeled_align

# Step 8 aligning files to be trimmed (muted) and then RMSD (Do for each pipeline) MAY NOT NEED TO DO IF WE DO 8.5
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled/RepBuilder/90/"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled/RepBuilder/80/"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled/TCRModel/90/"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled/TCRModel/80/"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled/LYRA/90/"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled/LYRA/80/"
# print('SeqA: Crystal Seq, SeqB: Modeled Seq')
# tool.unmute_all('/Users/austinseamann/PycharmProjects/TCRDOCK/Testing_align/')
# error_list = []
# for filename in os.listdir(directory):
#     if filename.endswith(".pdb"):
#         tool.set_file_name(directory[45:] + filename)  # Sets file name with removing first part of directory
#         tool.semi_align_pdb(filename[:4])
#         tool.set_file_name('Testing_align/' + filename[:4] + '_tcr.pdb')
# print(error_list)
# print('Error #: ' + str(len(error_list)))

# Step 8.25: Update count on modeled structures NEED TO STILL DO FOR REST OF SEQS.
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled_align/RepBuilder/90/"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled_align/RepBuilder/80/"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled_align/TCRModel/90/"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled_align/TCRModel/80/"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled_align/LYRA/90/"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled_align/LYRA/80/"
# for pdb in os.listdir(directory):
#     if pdb.endswith('.pdb'):
#         tool.set_file_name(directory[45:] + pdb)
#         print(tool.get_file_name())
#         tool.recount_tcr()

# Step 8.5 aligning files to be trimmed (muted) and using longest common sequences between all three modeling pipeline
# tool.lcs_pdbs(80)  # Switch between percentage

# Step 8.75 muting with lcs USE THIS ONE!!!
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled_align/RepBuilder/90/"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled_align/RepBuilder/80/"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled_align/TCRModel/90/"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled_align/TCRModel/80/"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled_align/LYRA/90/"
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled_align/LYRA/80/"
# tool.unmute_all(directory)
# tool.unmute_all("/Users/austinseamann/PycharmProjects/TCRDOCK/Testing_align/")
# error_list = []
# tool.lcs_semi_align('lcs_result_80.txt', 80, 'LYRA')

# Step 9
# RMSD time! Go over to Ubuntu when you finish modeling.

# Step 10 Make CDR files (Have to unmute all in Modeld_align)
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Testing_full"
# cdr_list = tool.make_cdr_ref()
# for each in sorted(os.listdir(directory)):
#     if each.endswith(".pdb"):
#         tool.set_file_name(directory + "/" + each)
#         tool.create_zone_entry()

# Check seqs
# directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/Testing_align"
# print('SeqA: Template Seq, SeqB: Modeled Seq')
# for filename in os.listdir(directory):
#     if filename.endswith(".pdb"):
#         tool.set_file_name('Testing_align/' + filename)
#         tool.check_seq('Modeled/RepBuilder/80')

# tool.set_file_name('Modeled/RepBuilder/90/1bd2.pdb')
# tool.clean_repbuilder()


# print(tool.prepare_seq('s2a'))
# print(tool.prepare_seq('s3a'))
# model_list = ['ALL_PDB/1bd2.pdb']
# for pdb in model_list:
#     tool.set_file_name(pdb)
#     print(tool.get_pdb_id())
#     print(tool.get_tcr_chains())
#     tool.spilt_tcr()

# tool.align_tcr_global('ALL_PDB/4g8g.pdb', 'ALL_PDB/3gsn.pdb')

# Use for finding best score between tcrs
# with open('PDB_NAMES.txt', 'r') as f:
#     name = ''
#     test = 'ALL_PDB/4g8g.pdb'
#     save = []
#     for line in f:
#         try:
#             name = 'ALL_PDB/' + line[:-1] + '.pdb'
#             score = tool.align_tcr_global(test, name)
#             save.append([score, name])
#         except:
#             print('Bad')
#     print(sorted(save))

# Use to confirm that each PDB only has two chains in file
# directory = '/Users/austinseamann/PycharmProjects/TCRDOCK/STCRDab_abTCR_clean'
# for filename in os.listdir(directory):
#     if filename.endswith('.pdb'):
#         tool.set_file_name('STCRDab_abTCR_clean/' + filename)
#         count = len(tool.get_chains())
#         if count != 2:
#             print("================")
#             print(filename)
#             print(tool.get_tcr_chains())
#             print("================")
#         else:
#             print(tool.get_file_name())
#             print(tool.get_tcr_chains())


# Use for creating the GENE_FAMILY Folder
# with open('PDB_NAMES.txt', 'r') as f:
#     name = ''
#     for line in f:
#         try:
#             name = "ALL_PDB/" + line[:-1] + ".pdb"
#             tool.file_name = name
#             print(name)
#             print(tool.get_tcr_chains())
#             tool.group_gene_family()
#         except:
#             print('Bad')

tool.set_file_name("/Users/austinseamann/PycharmProjects/TCRDOCK/Fran_Breakdown/Models/ES556M1BM-09.pdb")
print(tool.get_tcr_chains())
print(tool.get_amino_acid_on_chain("A"))
print(tool.get_amino_acid_on_chain("B"))
