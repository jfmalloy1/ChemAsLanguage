from rdkit import Chem
import pandas as pd
from tqdm import tqdm
import os
import numpy as np
from scipy import stats
import pickle
import random
from multiprocessing import Pool
import time
import json

""" Read a file containing fragments & occurrence counts
    Input: filepath
    Output: pandas dataframe with two columns, ["Frags", "Occurances"]
"""
def read_occurrences(fp):
    return pd.read_csv(fp)

""" Deletes all fragments found above the 1-quantile (still trying to figure out what this means)
    Input: dataframe (assumed to have an "Occurances" column), quantile q
    Output: filtered dataframe
"""
def delete_promiscuous_frags(df, q):
    return df[df["Occurances"] <= df["Occurances"].quantile(q)]

""" Link compound classes with specific compounds
    Input: none (assumes br08001.json as a file in working directory)
    Output: Dictionary of KEGG cpd ids associated with cpd class
"""
def find_cpd_classes():
    fp = "br08001.json"
    with open(fp) as json_file:
        br_json = json.load(json_file)
    #Create a dictionary with each compound having the appropriate label
    classes = []
    cpd_class = {}
    for key1 in br_json["children"]:
      classes.append(key1["name"])
      for key2 in key1["children"]:
        for key3 in key2["children"]:
          for cpd in key3["children"]:
            cpd_class[cpd["name"][:6] + ".mol"] = key1["name"]

    return cpd_class

""" Find compounds associated with a particular compound class (e.g., Lipids, Organic Acids, etc...)
    Input: filepath to specific json file containing KEGG information
    Output: smiles string associated with a particular cpd class
"""
def get_cpd_class():
    print(find_cpd_classes())

""" Calculate the percentage of all molecules containing a specific fragment
    Input: fragment, list of all compounds
    Output: percentage of compounds where a given fragment appears
"""
def coverage(frag, cpds):
    pass

def main():
    ### Core Fragments within Lipids ###
    #Test of core fragments
    lipid_df = read_occurrences("Biology/Data/Tests/Lipids_sampleOccurances/0_occurances.csv")
    print(lipid_df.head())

    #Delete "promiscuous" MCS fragments - all molecules aboves the 1st quantile (defined by F&M 2021)
    print(stats.describe(lipid_df["Occurances"]))
    print(lipid_df["Occurances"].quantile(0.25)) #Open Q - is 1st quantile 25%?
    print("Num of non-promiscuous frags:", len([x for x in lipid_df["Occurances"] if x <= lipid_df["Occurances"].quantile(0.25)]))
    filtered_lipid_df = delete_promiscuous_frags(lipid_df, 0.25)

    #Select candidate fragments - 3 step process
    #Step 1: coverage (% of molecules containing each MCS) - using OG molecule list
    #TODO: finish this! (use find_domains.py as a guide)
    get_cpd_class()

    #Step 2: homogeneity (avg pairwise Tanimoto similarity between all molecules sharing a MCS frag)

    #Step 3: inclusion (% of all other MCS frags which are substructures of a given MCS frag)

    #Step 4: MCScore = coverage * homogeneity * inclusion

    #Step 5: filter MCScore (keep only frags >= 70-quantile (70th percentile?))


if __name__ == "__main__":
    main()
