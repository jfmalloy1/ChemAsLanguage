import pandas as pd
from rdkit import Chem
import numpy as np
import re
import itertools

""" Return a dataframe of all KEGG data """
""" Input: the filepath (fp) to a csv file containing canonical smiles strings """
""" Output: dataframe of KEGG data, including canoncial smiles strings """
def get_KEGG_smiles(fp):
    df = pd.read_csv(fp)
    return df

""" Return a list of all chemical fragments in SMARTS form """
""" Input: filepath (fp) to txt file """
""" Output: list of SMARTS fragments """
def get_chem_fragments(fp):
    fragments = []
    with open(fp, "r") as f:
        for line in f:
            fragments.append(line.strip())

    return fragments

""" Determine the presence of a chemical fragment within a smiles string """
""" Input: smiles string of a chemical, fragment of a molecule in SMARTS format """
""" Output: T/F (presence/absence) """
#def find_single_presence(smiles_mol, f1_mol, f2_mol):
def find_single_presence(f1_occur, f2_occur):
    return (f1_occur and f2_occur)

def main():
    #Read in KEGG data
    df = get_KEGG_smiles("kegg_data.csv")
    #Read in chem fragments
    fragments = get_chem_fragments("frags.txt")

    #Co-occurance dataframe
    co_df = pd.DataFrame(0, index=fragments, columns=fragments)

    #Determine Mol structure of all KEGG (speed purposes)
    KEGG_mol_dict = {}
    for item, row in df.iterrows():
        KEGG_mol_dict[row["MOL file"]] = Chem.MolFromSmiles(row["Original SMILES"])

    #Loop through all fragments to find co-occurances
    count = 0
    occurance_dict = {}
    for f in fragments:
        count += 1
        print(count, "of", len(fragments))
        cpd_dict = {}
        for cpd in KEGG_mol_dict:
            try:
                cpd_dict[cpd] = KEGG_mol_dict[cpd].HasSubstructMatch(Chem.MolFromSmarts(f))
            except:
                cpd_dict[cpd] = False
        occurance_dict[f] = cpd_dict

    count = 0
    percent = -1
    #Get all pairs
    frag_combos = list(itertools.combinations(fragments, 2))
    l = len(frag_combos)
    for pair in frag_combos:
        #Output progress
        count += 1
        p = int(count / float(l) * 100)
        if p > percent:
            percent = p
            print(percent, "completed!")

        #Build co-occurance matrix
        for cpd in KEGG_mol_dict:
            if find_single_presence(occurance_dict[pair[0]][cpd], occurance_dict[pair[1]][cpd]):
                #If true, add one to dataframe count
                co_df.loc[pair[0]][pair[1]] += 1
                co_df.loc[pair[1]][pair[0]] += 1

    print(co_df)
    co_df.to_csv("frags_coocurrance.csv")

if __name__ == "__main__":
    main()
