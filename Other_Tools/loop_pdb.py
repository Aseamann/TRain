import chimera
import os
import Midas

directory = "/Users/austinseamann/PycharmProjects/TCRDOCK/"
folder = 'STCRDab_10_10_abTCR/'
location = directory + folder

for filename in sorted(os.listdir(location)):
    if filename.endswith('.pdb'):
        chimera.openModels.open(location + filename)
        Midas.focus('#')
        i = raw_input('Press ENTER for next PDB\nPress c and ENTER to exit')
        if i == 'c':
            break
        while i != '' or i == 'c':
            if i == 'c':
                break
            i = raw_input('Press ENTER for next PDB')
        chimera.closeSession()
print('No more PDBs')

