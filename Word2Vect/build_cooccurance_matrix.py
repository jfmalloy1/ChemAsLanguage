import pandas as pd
from rdkit import Chem
import numpy as np
import re

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
def find_single_presence(smiles_mol, f1_mol, f2_mol):
    try:
        return (smiles_mol.HasSubstructMatch(f1_mol) and smiles_mol.HasSubstructMatch(f2_mol))
    except:
        return False

def main():
    #Read in KEGG data
    df = get_KEGG_smiles("kegg_data.csv")
    #Read in chem fragments
    fragments = get_chem_fragments("functional_groups_smarts.txt")

    #Co-occurance dataframe
    co_df = pd.DataFrame(0, index=fragments, columns=fragments)

    #Determine Mol structure of all KEGG (speed purposes)
    KEGG_mol_dict = {}
    for item, row in df.iterrows():
        KEGG_mol_dict[row["MOL file"]] = Chem.MolFromSmiles(row["Original SMILES"])

    #Loop through all fragments to find co-occurances
    count = 0
    for f1 in fragments:
        count += 1
        print(f1, count, "of 314")
        #Mol structure for initial comparison fragment
        frag1_mol = Chem.MolFromSmarts(f1)
        #Compare each fragment with each other
        for f2 in fragments:
            #2nd fragment Mol structure
            frag2_mol = Chem.MolFromSmarts(f2)
            #Loop through all KEGG cpds
            for cpd in KEGG_mol_dict:
                if find_single_presence(KEGG_mol_dict[cpd], frag1_mol, frag2_mol):
                    #If true, add one to dataframe count
                    co_df.loc[f1][f2] += 1



    print(co_df)
    co_df.to_csv("functional_groups_coocurrance.csv")

if __name__ == "__main__":
    main()
