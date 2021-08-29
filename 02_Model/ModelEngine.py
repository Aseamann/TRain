# This file is a part of the TRain program
# Author: Austin Seamann & Dario Ghersi
# Version: 0.01
# Last Updated: Aug 23rd, 2021
import argparse
import subprocess
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time
import os


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


def lyra_run(fasta, start_pdb):
    error_list = []
    with open(fasta, 'r') as file:
        pdb_list = {}
        current_pdb = ""
        alpha = False
        beta = False
        blacklist = False
        start = False
        for line in file:
            if line[:-1] == start_pdb or start_pdb == "...":
                start = True
            if start:
                if alpha:
                    pdb_list[current_pdb]['ALPHA'] = line[:-1]
                    alpha = False
                elif beta:
                    pdb_list[current_pdb]['BETA'] = line[:-1]
                    beta = False
                elif blacklist:
                    pdb_list[current_pdb]["Blacklist"] = line[:-1]
                    blacklist = False
                elif line.startswith("ALPHA:"):
                    alpha = True
                elif line.startswith("BETA:"):
                    beta = True
                elif line.startswith("Black"):
                    blacklist = True
                elif not alpha and not beta and not blacklist:
                    if line != "\n":
                        pdb_list[line[:-1]] = {}
                        current_pdb = line[:-1]
        print(pdb_list)
    for entry in pdb_list:
        driver.get("http://www.cbs.dtu.dk/services/LYRA/index.php")
        try:
            element = WebDriverWait(driver, 10).until(
                EC.presence_of_element_located((By.ID, "chain1"))
            )

            pdb_id = entry + '.pdb'
            chain1 = driver.find_element_by_id("chain1")
            chain1.send_keys(pdb_list[entry]['ALPHA'])

            chain2 = driver.find_element_by_id("chain2")
            chain2.send_keys(pdb_list[entry]['BETA'])

            blacklist = driver.find_element_by_id("blacklistPDBs")
            blacklist.send_keys(pdb_list[entry]['Blacklist'])

            time.sleep(1)
            submit = driver.find_elements_by_xpath("//input[@value='Submit']")[0]
            submit.click()

            try:
                element = WebDriverWait(driver, 100).until(
                    EC.presence_of_element_located((By.ID, "jsmol"))
                )
                time.sleep(3)
                download = driver.find_elements_by_xpath("/html/body/div[1]/div[2]/div[2]/nav/div[2]/ul[2]/div/a")[0]
                download.click()
                time.sleep(5)
                for pdb in os.listdir(download_path):
                    if pdb.startswith("lyra"):
                        os.rename(download_path + "/" + pdb, download_path + "/" + pdb_id)
                print(pdb_id)
            except:
                error_list.append(pdb_id)
                print("Error: Level 1: " + pdb_id)

        except:
            print("Error: Level 0")


def tcrmodel_run(fasta, start_pdb):
    error_list = []
    with open(fasta, 'r') as file:
        pdb_list = {}
        current_pdb = ""
        alpha = False
        beta = False
        blacklist = False
        start = False
        for line in file:
            if line[:-1] == start_pdb or start_pdb == "...":
                start = True
            if start:
                if alpha:
                    pdb_list[current_pdb]['ALPHA'] = line[:-1]
                    alpha = False
                elif beta:
                    pdb_list[current_pdb]['BETA'] = line[:-1]
                    beta = False
                elif blacklist:
                    pdb_list[current_pdb]["Blacklist"] = line[:-1]
                    blacklist = False
                elif line.startswith("ALPHA:"):
                    alpha = True
                elif line.startswith("BETA:"):
                    beta = True
                elif line.startswith("Black"):
                    blacklist = True
                elif not alpha and not beta and not blacklist:
                    if line != "\n":
                        pdb_list[line[:-1]] = {}
                        current_pdb = line[:-1]
        print(pdb_list)

    for entry in pdb_list:
        driver.get("https://tcrmodel.ibbr.umd.edu/index")
        try:
            element = WebDriverWait(driver, 10).until(
                EC.presence_of_element_located((By.ID, "formid"))
            )

            pdb_id = entry + '.pdb'
            chain1 = driver.find_element_by_id("achain")
            chain1.send_keys(pdb_list[entry]['ALPHA'])

            chain2 = driver.find_element_by_id("bchain")
            chain2.send_keys(pdb_list[entry]['BETA'])

            advanced = driver.find_elements_by_xpath("/html/body/div/article/div[2]/form/div[5]/u")[0]
            advanced.click()

            time.sleep(1)
            blacklist = driver.find_elements_by_name("pdbblacklist")[0]
            blacklist.clear()
            blacklist.send_keys(pdb_list[entry]['Blacklist'])

            time.sleep(1)
            submit = driver.find_elements_by_xpath("/html/body/div/article/div[2]/form/div[6]/input[3]")[0]
            submit.click()

            try:
                flag = False
                time.sleep(2)
                while loading("/html/body/div/article/h4"):
                    time.sleep(2)
                if error_path("/html/body/div/article/h4") != True:
                    element = WebDriverWait(driver, 20).until(
                        EC.presence_of_element_located((By.ID, "download-wrap"))
                    )
                    time.sleep(3)
                    download = driver.find_elements_by_xpath("/html/body/div/article/div[1]/div[3]/h4[3]/a")[0]
                    download.click()
                    time.sleep(5)
                    for pdb in os.listdir(download_path):
                        if pdb.startswith("tcrmodel"):
                            os.rename(download_path + "/" + pdb, download_path + "/" + pdb_id)
                    print(pdb_id)
                else:
                    error_list.append(pdb_id)
                    print("Error: Level 2: " + pdb_id)
            except:
                error_list.append(pdb_id)
                print("Error: Level 1: " + pdb_id)
        except:
            print("Error: Level 0")


def loading(xpath):
    try:
        text = driver.find_elements_by_xpath(xpath)[0].text
        if str(text).__contains__('submitted'):
            return True
        else:
            return False
    except:
        return False


def error_path(xpath):
    try:
        text = driver.find_elements_by_xpath(xpath)[0].text
        if str(text).__contains__('Failed'):
            return True
        else:
            return False
    except:
        return False


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="Fasta file containing TCR sequeces for modeling.", type=str)
    parser.add_argument("-a", "--alpha", help="Fasta file containing each alpha chain sequence of TCR. Must match "\
                        "Beta header.", type=str)
    parser.add_argument("-b", "--beta", help="Fasta file containing each beta chain sequence of TCR. Must match "\
                        "Alpha header.", type=str)
    parser.add_argument("-l", "--lyra", help="Automate LYRA submission", default=False, action="store_true")
    parser.add_argument("-t", "--tcrmodel", help="Automate TCRmodel submission", default=False,
                        action="store_true")
    parser.add_argument("-s", "--start", help="Optional starting position in fasta file. Provide header.", type=str,
                        default="...")
    parser.add_argument("-r", "--results", help="Download directory folder name", default="Models", type=str)
    return parser.parse_args()


def main():
    args = parse_args()
    modeling_input(fasta, args)
    chromeOptions = webdriver.ChromeOptions()
    os.mkdir(args.results)
    download_path = os.getcwd() + "/" + args.results
    prefs = {"download.default_directory": download_path}
    chromeOptions.add_experimental_option("prefs", prefs)
    chromedriver = "/Applications/chromedriver"
    driver = webdriver.Chrome(executable_path=chromedriver, options=chromeOptions)
    if args.lyra:
        lyra_run(args.fasta, args.start)
    if args.tcrmodel:
        tcrmodel_run()
    time.sleep(30)
    driver.quit()


if __name__ == '__main__':
    main()
