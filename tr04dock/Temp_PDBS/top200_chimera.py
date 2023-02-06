import chimera
import os
import Midas

directory = os.getcwd() + "/"
path = directory + "top200/"

flag = True

chimera_colors = ["red", "yellow", "green", "cyan", "light sea green", "blue", "orange",
                  "cornflower blue", "purple", "forest green", "hot pink", "magenta",
                  "orange red" "gray", "dim gray"]
pdbs = []
cluster = []

with open(path + "top200.csv", "r") as f:
    for line in f:
        print(line)
        split = line.split("\t")
        print(split)
        if split[0] != "PDB":
            pdbs.append(split[0] + ".pdb")
            cluster.append(int(split[1]))

num_colors = max(cluster)

for i in range(len(pdbs)):
    chimera.openModels.open(path + "/" + pdbs[i])
    if flag:
        Midas.focus('#')
        flag = False
    print(chimera_colors[cluster[i]])
    chimera_cmd = "color " + chimera_colors[cluster[i]] + ",r,l,s #" + str(i)
    print(chimera_cmd)
    chimera.runCommand(chimera_cmd)
print('No more PDBs')

