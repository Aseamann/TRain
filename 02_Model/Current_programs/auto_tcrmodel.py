from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time
import os

chromeOptions = webdriver.ChromeOptions()
download_path = "/Users/austinseamann/PycharmProjects/TCRDOCK/Modeled/TCRModel/80"
prefs = {"download.default_directory" : download_path}
chromeOptions.add_experimental_option("prefs",prefs)
chromedriver = "/Applications/chromedriver"
driver = webdriver.Chrome(executable_path=chromedriver, options=chromeOptions)

starting_pdb = '4grm'
error_list = []


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

print("Complete")
print(error_list)
time.sleep(30)
driver.quit()

