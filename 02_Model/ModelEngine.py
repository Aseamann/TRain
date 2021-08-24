# This file is a part of the TRain program
# Author: Austin Seamann & Dario Ghersi
# Version: 0.01
# Last Updated: Aug 23rd, 2021
import argparse
import subprocess


# TODO: write dependencies of cd-hit
def check_redundant(fasta, percent):
    subprocess.run(["cd-hit", "-i", fasta, "-o", "clusters.txt", "-c", percent, "-n", 5, "-M", 8000, "-d", 0])


def modeling_input(fasta, service):
    if service.rep_builder:
        with open("alpha.fasta", "w+") as f1:
            with open("beta.fasta", "w+") as f2:
                print("help")
    if service.lyra:
        print("help")
    if service.tcrmodel:
        print("help")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="Fasta file containing TCR sequeces for modeling.", type=str)
    parser.add_argument("-r", "--rep_builder", help="(Output) Repertoire Builder input file format", default=False,
                        action="store_true")
    parser.add_argument("-l", "--lyra", help="(Output) Lyra input file format", default=False, action="store_true")
    parser.add_argument("-t", "--tcrmodel", help="(Output) TCRmodel input file format", default=False,
                        action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    modeling_input(fasta, args)


if __name__ == '__main__':
    main()