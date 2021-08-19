import os
import subprocess

directory = os.getcwd()
for file in os.listdir(directory):
    content = subprocess.run("tr -s ' ' < " + file + " | tr ' ' ','", shell=True, stdout=subprocess.PIPE)
    with open(file[:-4] + ".csv", "w") as f:
        for line in content.stdout.decode('utf-8').split('\n'):
            f.write(line + '\n')
