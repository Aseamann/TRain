# This file is a part of the TRain program
# Author: Austin Seamann & Dario Ghersi
# Version: 0.1
# Last Updated: January 12th, 2022
import argparse
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time
import os


# Read in alpha and beta fasta files and generate tcr_seqs dictionary
def read_fastas(alpha_in, beta_in, start_pdb):
    # Dictionary created (e.g. tcr_seqs[CHAIN][HEADER] = seq)
    list_files = {"ALPHA": alpha_in, "BETA": beta_in}
    tcr_seqs = {"ALPHA": {}, "BETA": {}}
    if start_pdb != "...":  # Starting further down in fasta file
        start = False
    else:
        start = True
    for chain in list_files:
        with open(list_files[chain], "r") as f1:
            header = ""
            for line in f1:
                if not start:  # Searching for first header
                    if line.split(">")[:-1].split(">")[1] == header:
                        start = True
                if start:  # Wont start until first header is found
                    if line[0] == ">":
                        header = line[:-1].split(">")[1]
                    elif header != "":  # if we've started our first header
                        if header in tcr_seqs[chain].keys():  # If string already started
                            tcr_seqs[chain][header] += line[:-1]
                        else:  # If first line of seq
                            tcr_seqs[chain][header] = line[:-1]
    return tcr_seqs


# Automation of LYRA website submission of alpha and beta chains to download to folder
def lyra_run(tcr_seqs, driver, download_path, scalar):
    error_list = []
    for header in tcr_seqs["ALPHA"]:
        driver.get("https://services.healthtech.dtu.dk/service.php?LYRA-1.0")
        num = 10 * scalar
        time.sleep(num)
        driver.switch_to.frame(driver.find_element_by_tag_name("iframe"))
        try:
            num = 10 * scalar
            WebDriverWait(driver, num).until(
                EC.presence_of_element_located((By.CSS_SELECTOR, '#chain1'))
            )
            num = 5 * scalar
            time.sleep(num)

            if "(" in header:  # Detect if there's a blacklist in header
                parts = header.split("((")
                header = parts[0].strip()
                blacklist = parts[1].strip("()").strip
                driver.find_element_by_xpath('//*[@id="blacklistPDBs"]').send_keys(blacklist)
            pdb_id = header + '.pdb'
            #chain1 = driver.find_element_by_id("chain1")
            chain1 = driver.find_elements_by_xpath('//*[@id="chain1"]')[0]
            chain1.send_keys(tcr_seqs["ALPHA"][header])

            #chain2 = driver.find_element_by_id("chain2")
            chain2 = driver.find_elements_by_xpath('//*[@id="chain2"]')[0]
            chain2.send_keys(tcr_seqs["BETA"][header])

            num = 1 * scalar
            time.sleep(num)
            submit = driver.find_elements_by_xpath('/html/body/div/div/div/form/div[4]/div/input')[0]
            submit.click()

            try:
                num = 200 * scalar
                WebDriverWait(driver, num).until(
                    EC.presence_of_element_located((By.ID, "jsmol"))
                )
                num = 3 * scalar
                time.sleep(num)
                download = driver.find_elements_by_xpath("/html/body/div[1]/div[2]/div[2]/nav/div[2]/ul[2]/div/a")[0]
                download.click()
                num = 5 * scalar
                time.sleep(num)
                for pdb in os.listdir(download_path):
                    if pdb.startswith("lyra"):
                        os.rename(download_path + "/" + pdb, download_path + "/" + pdb_id)
                print(pdb_id)
            except:
                error_list.append(pdb_id)
                print("Error: Level 1: " + pdb_id)

        except:
            print("Error: Level 0")
    time.sleep(30)
    driver.quit()


def tcrmodel_run(tcr_seqs, driver, download_path, scalar):
    error_list = []
    for header in tcr_seqs["ALPHA"]:
        driver.get("https://tcrmodel.ibbr.umd.edu/index")
        try:
            num = 10 * scalar
            WebDriverWait(driver, num).until(
                EC.presence_of_element_located((By.ID, "formid"))
            )

            pdb_id = header + '.pdb'
            chain1 = driver.find_element_by_id("achain")
            chain1.send_keys(tcr_seqs["ALPHA"][header])

            chain2 = driver.find_element_by_id("bchain")
            chain2.send_keys(tcr_seqs["BETA"][header])

            advanced = driver.find_elements_by_xpath("/html/body/div/article/div[2]/form/div[5]/u")[0]
            advanced.click()

            num = 1 * scalar
            time.sleep(num)
            submit = driver.find_elements_by_xpath("/html/body/div/article/div[2]/form/div[6]/input[3]")[0]
            submit.click()

            try:
                num = 2 * scalar
                time.sleep(num)
                while loading("/html/body/div/article/h4", driver):
                    time.sleep(num)
                if not error_path("/html/body/div/article/h4", driver):
                    num = 200 * scalar
                    element = WebDriverWait(driver, num).until(
                        EC.presence_of_element_located((By.ID, "download-wrap"))
                    )
                    num = 3 * scalar
                    time.sleep(num)
                    download = driver.find_elements_by_xpath("/html/body/div/article/div[1]/div[3]/h4[3]/a")[0]
                    download.click()
                    num = 5 * scalar
                    time.sleep(num)
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
    time.sleep(30)
    driver.quit()


# Automation of TCRBuilder website submission of alpha and beta chains to download to folder
def tcrbuilder(tcr_seqs, driver, download_path, scalar):
    error_list = []
    for header in tcr_seqs["ALPHA"]:
        driver.get("http://opig.stats.ox.ac.uk/webapps/stcrpred/TCRBuilder/")
        try:
            num = 10 * scalar
            WebDriverWait(driver, num).until(
                EC.presence_of_element_located((By.ID, "collapseTwo"))
            )

            if "(" in header:  # Detect if there's a blacklist in header
                parts = header.split("((")
                header = parts[0].strip()
                blacklist = parts[1].strip("()").strip
                driver.find_element_by_css_selector('#jobsubmission > div:nth-child(9) > textarea').send_keys(blacklist)
            pdb_id = header + '.pdb'
            chain1 = driver.find_element_by_css_selector('#jobsubmission > div:nth-child(3) > textarea')
            chain1.send_keys(tcr_seqs["ALPHA"][header])

            chain2 = driver.find_element_by_css_selector('#jobsubmission > div:nth-child(5) > textarea')
            chain2.send_keys(tcr_seqs["BETA"][header])

            num = 1 * scalar
            time.sleep(num)
            submit = driver.find_element_by_css_selector('#jobsubmission > button')
            submit.click()

            try:
                num = 1800 * scalar
                WebDriverWait(driver, num).until(
                    EC.text_to_be_present_in_element((By.CSS_SELECTOR,
                                                      '#jobstatus > table > tbody > tr:nth-child(3) > td:nth-child(2) > font'), "finished!")
                )
                num = 3 * scalar
                time.sleep(num)
                download = driver.find_element_by_css_selector(
                    '#downloads > table > tbody > tr:nth-child(4) > td:nth-child(2) > a'
                )
                download.send_keys(download_path)
                num = 5 * scalar
                time.sleep(num)
                temp_dir = os.getcwd()
                download = driver.find_elements_by_xpath('//*[@id="downloads"]/table/tbody/tr[2]/td[2]/a')[0]
                download.click()
                time.sleep(num)
                for pdb in os.listdir(download_path):
                    if pdb.startswith("fread"):
                        os.rename(download_path + "/" + pdb, download_path + "/" + pdb_id)
                    if pdb.endswith(".log"):
                        os.remove(pdb)
                        error_list.append(pdb_id)
                print(pdb_id)
            except:
                error_list.append(pdb_id)
                print("Error: Level 1: " + pdb_id)

        except:
            print("Error: Level 0")
    time.sleep(30)
    driver.quit()


def loading(xpath, driver):
    try:
        text = driver.find_elements_by_xpath(xpath)[0].text
        if str(text).__contains__('submitted'):
            return True
        else:
            return False
    except:
        return False


def error_path(xpath, driver):
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
    parser.add_argument("-a", "--alpha", help="Fasta file containing each alpha chain sequence of TCR. Must match "\
                        "Beta header.", type=str)
    parser.add_argument("-b", "--beta", help="Fasta file containing each beta chain sequence of TCR. Must match "\
                        "Alpha header.", type=str)
    parser.add_argument("-l", "--lyra", help="Automate LYRA submission", default=False, action="store_true")
    parser.add_argument("-m", "--tcrmodel", help="Automate TCRmodel submission", default=False,
                        action="store_true")
    parser.add_argument("-t", "--tcrbuilder", help="Automate TCRBuilder submission", default=False, action="store_true")
    parser.add_argument("-s", "--start", help="Optional starting position in fasta file. Provide header.", type=str,
                        default="...")
    parser.add_argument("-d", "--download", help="Download directory folder name", default="Models", type=str)
    parser.add_argument("-g", "--long", help="Scalar for wait times for elements to appear and page to load", default=1.0,
                        type=float)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.download not in os.listdir():
        os.mkdir(args.download)
    with open("config.ini", "r") as f1:  # Grab webdriver location
        for line in f1:
            if line[:10] == "driver_loc":
                chromedriver = line[:-1].split("=")[1][1:-1]
    chrome_options = webdriver.ChromeOptions()
    download_path = os.getcwd() + "/" + args.download
    prefs = {"download.default_directory": download_path}
    chrome_options.add_experimental_option("prefs", prefs)
    driver = webdriver.Chrome(executable_path=chromedriver, options=chrome_options)
    if args.lyra:
        # Create tcr_seqs dictionary containing alpha and beta chains
        tcr_seqs = read_fastas(args.alpha, args.beta, args.start)
        lyra_run(tcr_seqs, driver, download_path, args.long)
    if args.tcrmodel:
        # Create tcr_seqs dictionary containing alpha and beta chains
        tcr_seqs = read_fastas(args.alpha, args.beta, args.start)
        tcrmodel_run(tcr_seqs, driver, download_path, args.alpha, args.beta, args.long)
    if args.tcrbuilder:
        # Create tcr_seqs dictionary containing alpha and beta chains
        tcr_seqs = read_fastas(args.alpha, args.beta, args.start)
        tcrbuilder(tcr_seqs, driver, download_path, args.long)


if __name__ == '__main__':
    main()
