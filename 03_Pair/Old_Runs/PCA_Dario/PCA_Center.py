import argparse
import statistics
import numpy as np
from numpy import linalg
from sklearn.decomposition import PCA
from scipy.spatial.transform import Rotation
from math import sqrt


class PCAcenter:
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

    # Returns a list of all chains in PDB file
    def get_chains(self):
        chains = []
        with open(self.file_name, 'r') as file:
            for line in file:
                if line[0:6] == 'ATOM  ':
                    if not chains.__contains__(line[21]):
                        chains.append(line[21])
        return chains

    # Collect the atoms from an inputted chain. Provides all values in PDB file
    def get_atoms_on_chain(self, chain):
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

    # Rebuilt lines for atoms information submitted
    def rebuild_atom_line(self, atoms):
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
        # Determinate of components matrix - Corrects if determinate is negative
        determinate = np.linalg.det(pca.components_)
        print("determinate")
        if determinate < 0:
            for position in pca.components_:
                position[0] = position[0] * -1
        transpose = np.transpose(pca.components_)
        ####################################################################################
        #   NEW    #
        ############
        # Always have N-terminus in positive coordinates (Fixes flips on x-axis)
        test_x = np.matmul(new_array[-1], transpose)
        if test_x[1] < 0:
            print("X-axis")
            r = Rotation.from_euler('x', 180, degrees=True)
            transpose = r.apply(transpose)
        # Always have first chain on left side | alpha on left side (Fixes flips on y-axis)
        # Determines by looking at last atom which should be on the Beta chain
        test_y = np.matmul(new_array[-1], transpose)
        print(test_y)
        if test_y[0] < 0:
            print("Y-axis")
            r = Rotation.from_euler('y', 180, degrees=True)
            transpose = r.apply(transpose)
        #####################################################################################
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
        # --- Bug testing ---Print initial center coord.
        print("Previous Center")
        print(self.get_center())
        # Print updated center coord.
        print("Center")
        self.set_file_name(new_name)
        print(self.get_center())


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="PDB file input", type=str)
    return parser.parse_args()


def main():
    args = parse_args()
    pdb = PCAcenter(args.pdb)
    pdb.center()


if __name__ == '__main__':
    main()
