from rdkit import Chem
import pandas as pd
from tqdm import tqdm
from scipy import stats
import pickle

""" Get smiles strings for a set of compounds
    Input: list of cpd ids for a given set of compounds
    Output: the smiles strings associated with those compounds
"""
def get_smiles(fp):
    df = pd.read_csv(fp)
    kegg_df = pd.read_csv("Biology/Data/KEGG_chiral_molweight_formula_labels.csv")

    #Isolate archea SMILES (assumption: def has a column named "compounds")
    return kegg_df[kegg_df["C"].isin(df["compounds"].tolist())]["S"].tolist()

""" Convert a smiles string to rdkit mol objects
    Input: list of smiles string
    Output: list of RDkit mol objects
"""
def convert_to_mol_smiles(smiles):
    return [Chem.MolFromSmiles(smi.strip()) for smi in smiles]

""" Convert a smarts string to rdkit mol objects
    Input: list of smart strings
    Output: list of RDkit mol objects
"""
def convert_to_mol_smarts(smarts):
    return [Chem.MolFromSmarts(smt.strip()) for smt in smarts]

""" Find fragments (& associated statistics) within each compounds
    Input: A single cpd, entire list of Fragments
    Output: Number of fragments found within a cpd
"""
def find_frags_in_cpds(cpd, frags):
    frag_count = 0
    for f in frags:
        if cpd.HasSubstructMatch(f):
            frag_count += 1
    return frag_count

def main():
    #Read in fragments (Archaea as a test)
    fp = "Biology/Data/Archaea/Archaea0_fullOccurances.csv"
    a_df = pd.read_csv(fp)
    a_frags = a_df["Frags"].tolist() #Smarts strings found in archaea
    a_frags = convert_to_mol_smarts(a_frags)
    print("Found", len(a_frags), "archaea fragments")

    #Archaea smiles
    fp = "Biology/Data/archaea_cpds.csv"
    a_cpds = get_smiles(fp) #Smiles strings of archaea cpds
    a_cpds = [c for c in a_cpds if str(c) != 'nan'] #Remove nans
    a_cpds = convert_to_mol_smiles(a_cpds)
    print("Found", len(a_cpds), "archea compounds")

    #Find number of fragments in each compounds
    frag_count = []
    for cpd in tqdm(a_cpds):
        frag_count.append(find_frags_in_cpds(cpd, a_frags))

    #Basic descriptive statistics
    print(stats.describe(frag_count))
    pickle.dump(frag_count, open("Biology/Data/fragment_counts.p", "wb"))


if __name__ == "__main__":
    main()
