from rdkit import Chem
from rdkit import DataStructs
import pandas as pd
from tqdm import tqdm
import numpy as np
from scipy import stats
import pickle
from multiprocessing import Pool
import json
from itertools import combinations

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
    Input: cpd class label, assumes a file br08001.json containing cpd class info and Biology/Data/kegg_data.csv containing all KEGG data
    Output: RDKit mol objects associated with a particular cpd class, as well as cpd ids
"""
def get_cpd_class_mols(class_label):
    cpd_classes = find_cpd_classes()     #Data about all KEGG classes
    kegg_df = pd.read_csv("Biology/Data/kegg_data.csv") #KEGG data
    cpd_ids = [k for k,v in cpd_classes.items() if v == class_label] #Finding cpd labels for a class
    cpd_smiles = kegg_df[kegg_df["MOL file"].isin(cpd_ids)]["Original SMILES"].tolist() #Converting labels to smiles
    cpd_mols = [Chem.MolFromSmiles(x) for x in cpd_smiles] #Converting smiles to rdkit mol objects
    return cpd_mols, cpd_ids

""" Calculates the Tanimoto similarity between all cpds with the same fragment
    Input: list of RDKit mol objects
    Output: similarity score
"""
def similarity(mols):
    fps = [Chem.RDKFingerprint(x) for x in mols]
    similarity_score = []
    for c in combinations(fps,2):
        c1, c2 = c
        similarity_score.append(DataStructs.FingerprintSimilarity(c1,c2))

    return np.mean(similarity_score)

""" Calculate the percentage of all molecules containing a specific fragment
    Input: fragment, list of all compounds
    Output: percentage of compounds where a given fragment appears
"""
def get_coverage_and_similarity(frag, class_mols):
    frag = Chem.MolFromSmarts(frag)
    count = 0
    mols_with_frag = []
    for cm in class_mols:
        if cm.HasSubstructMatch(frag):
            count += 1
            mols_with_frag.append(cm)

    return count / len(class_mols), similarity(mols_with_frag)

""" Calculate the percentage of frags which are substructures of a given frag
    Input: fragment list
    Output: inclusion list
"""
def get_inclusion(frags):
    inclusion_list = []
    for frag in frags:
        inclusion_count = 0
        m = Chem.MolFromSmarts(frag)
        for f in frags:
            if f != frag and m.HasSubstructMatch(Chem.MolFromSmarts(f)):
                inclusion_count += 1
        inclusion_list.append(inclusion_count / len(frags))

    return inclusion_list

def main():
    ### Core Fragments within Compound Classes ###
    #Test of core fragments
    df = read_occurrences("Biology/Data/Tests/Carbs_sampleOccurances/0_occurances.csv")

    #Delete "promiscuous" MCS fragments - all molecules aboves the 1st quantile (defined by F&M 2021)
    print("Full occurances:", stats.describe(df["Occurances"]))
    print(df["Occurances"].quantile(0.25)) #Open Q - is 1st quantile 25%?
    print("Num of non-promiscuous frags:", len([x for x in df["Occurances"] if x <= df["Occurances"].quantile(0.25)]))
    print()
    filtered_df = delete_promiscuous_frags(df, 0.25)

    #Select candidate fragments - 3 step process
    #Step 1: coverage (% of molecules containing each MCS) - using OG molecule list
    #Step 2: homogeneity (avg pairwise Tanimoto similarity between all molecules sharing a MCS frag)

    class_mols, class_ids = get_cpd_class_mols("Carbohydrates")
    coverage_list = []
    similarity_list = []
    for frag in filtered_df["Frags"]:
        coverage_score, similarity_score = get_coverage_and_similarity(frag, class_mols)
        coverage_list.append(coverage_score)
        similarity_list.append(similarity_score)
    print("Coverage %:", stats.describe(coverage_list))
    print("Similarity:", stats.describe(similarity_list))

    #Step 3: inclusion (% of all other MCS frags which are substructures of a given MCS frag)
    inclusion_list = get_inclusion(filtered_df["Frags"])
    print("Inclusion %:", stats.describe(inclusion_list))

    #Step 4: MCScore = coverage * homogeneity * inclusion
    MCScore = np.prod(np.vstack([coverage_list, similarity_list, inclusion_list]), axis=0)
    print("MCSscore:", stats.describe(MCScore))
    print()

    #Step 5: filter MCScore (keep only frags >= 70-quantile (70th percentile?))
    cutoff = np.percentile(MCScore, 70)
    frag_list = filtered_df["Frags"].tolist()
    core_frags = []
    for i in range(len(frag_list)):
        if MCScore[i] >= cutoff:
            core_frags.append(frag_list[i])
    print("Number of core fragments:", len(core_frags))
    core_mols = [Chem.MolFromSmarts(x) for x in core_frags]
    print()

    ### Test core fragments over all KEGG ###
    #Create mol objects for each KEGG cpd
    kegg_df = pd.read_csv("Biology/Data/kegg_data.csv") #KEGG data
    print("Converting KEGG to RDKit mol objects...")
    kegg_df["Mols"] = [Chem.MolFromSmiles(x) for x in kegg_df["Original SMILES"]]
    print()

    #find all cpds where the core fragments appear
    print("Checking for core compounds...")
    core_cpds = []
    for index, row in tqdm(kegg_df.iterrows(), total=len(kegg_df)):
        for cmol in core_mols:
            try:
                if row["Mols"].HasSubstructMatch(cmol):
                    core_cpds.append(row["MOL file"])
            except:
                pass
    print()

    core_cpds = list(set(core_cpds))
    print("Number of core compounds:", len(core_cpds))

    print("Number of class compounds:", len(class_ids))

    overlapping_cpds = list(set(core_cpds).intersection(class_ids))
    print("Overlapping compounds:", len(overlapping_cpds))
    print("Percentage found in class cpds:", len([x for x in overlapping_cpds if x in class_ids]) / len(class_ids))




if __name__ == "__main__":
    main()
