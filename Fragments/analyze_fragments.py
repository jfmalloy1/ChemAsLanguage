from rdkit import Chem
import pandas as pd
from tqdm import tqdm
from scipy import stats
import pickle
import networkx as nx
from itertools import combinations
import matplotlib.pyplot as plt

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
    frags_in_cpd = []
    for f in frags:
        if cpd.HasSubstructMatch(f):
            frag_count += 1
            frags_in_cpd.append(Chem.MolToSmiles(f))
    return frag_count, frags_in_cpd

""" Create a networkx graph of fragments. Fragments are connected if they are within the same compound
    Input: List of lists, each sublist are fragments within a single archaea compound
    Output: A networkx graph of fragments
"""
def create_graph(frags):
    G = nx.Graph()
    #TEST: create a network of fragments from the 0th compound
    for cpd in tqdm(frags):
        for f in cpd:
            G.add_node(f)
        G.add_edges_from(combinations(cpd,2))
    print(nx.info(G))
    return G
    # #Try drawing it
    # plt.subplot(121)
    # nx.draw(G, with_labels=True, font_weight="bold")
    # plt.show()

def main():
    ### Number of fragments in each compound ###
    # #Read in fragments (Archaea as a test)
    # fp = "Biology/Data/Archaea/Archaea0_fullOccurances.csv"
    # a_df = pd.read_csv(fp)
    # a_frags = a_df["Frags"].tolist() #Smarts strings found in archaea
    # a_frags = convert_to_mol_smarts(a_frags)
    # print("Found", len(a_frags), "archaea fragments")
    #
    # #Archaea smiles
    # fp = "Biology/Data/archaea_cpds.csv"
    # a_cpds = get_smiles(fp) #Smiles strings of archaea cpds
    # a_cpds = [c for c in a_cpds if str(c) != 'nan'] #Remove nans
    # a_cpds = convert_to_mol_smiles(a_cpds)
    # print("Found", len(a_cpds), "archea compounds")
    #
    # #Find number of fragments in each compounds
    # frag_count = []
    # frags_in_all_cpds = []
    # for cpd in tqdm(a_cpds):
    #     count, frags_in_cpd = find_frags_in_cpds(cpd, a_frags)
    #     frag_count.append(count)
    #     frags_in_all_cpds.append(frags_in_cpd)
    #
    # #Basic descriptive statistics
    # print(stats.describe(frag_count))
    # pickle.dump(frag_count, open("Biology/Data/archaea_fragment_counts.p", "wb"))
    # pickle.dump(frags_in_all_cpds, open("Biology/Data/archaea_listOfAllFragsPerCpd.p", "wb"))

    ### Network of Compounds ###
    frags_in_all_cpds = pickle.load(open("Biology/Data/archaea_listOfAllFragsPerCpd.p", "rb"))

    pickle.dump(create_graph(frags_in_all_cpds), open("Biology/Data/archea_frag_network.p", "wb"))



if __name__ == "__main__":
    main()
