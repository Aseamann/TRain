import chimera
import os
import Midas

directory = "/Users/austinseamann/Library/Mobile Documents/com~apple~CloudDocs/School/BioI-Lab/TRain/Temp_Sample_Models/RepBuilder/"
folder = 'Results/'
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

