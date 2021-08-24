from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time
import os

chromeOptions = webdriver.ChromeOptions()
download_path = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled/LYRA/80"
prefs = {"download.default_directory" : download_path}
chromeOptions.add_experimental_option("prefs",prefs)
chromedriver = "/Applications/chromedriver"
driver = webdriver.Chrome(executable_path=chromedriver, options=chromeOptions)

starting_pdb = '1bd2'
error_list = []

with open("Modeled/all_80.txt", 'r') as file:
    pdb_list = {}
    current_pdb = ""
    alpha = False
    beta = False
    blacklist = False
    start = False
    for line in file:
        if line[:-1] == starting_pdb:
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

print("Complete")
print(error_list)
time.sleep(30)
driver.quit()

