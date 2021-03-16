from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
from tqdm import tqdm
from scipy import stats
import pickle
import networkx as nx
from itertools import combinations
import matplotlib.pyplot as plt
import os
import itertools

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

""" Get a list of mol objects from assembly fragments
    Input: File path (fp) containing assembly fragments
    Output: A list of RDKit mol objects
"""
def get_assembly_frags(fp):
    lines = open(fp).readlines()
    #Find beginning of inchi (next row after '$' character)
    count = 0
    for l in lines:
        if "$" in l:
            break
        count += 1

    #Turn all inchi fragments into mol objects
    return [Chem.MolFromInchi(i.strip()) for i in lines[count:] if "InChI" in i] #Note: Header takes 37 lines in this particular example

""" Convert a list of RDKit mol objects into (canonical) smiles
    Input: RDKit mol object list
    Output: list of canonical smiles
"""
def get_canonical_smiles(mols):
    smiles = list(set([Chem.MolToSmiles(m) for m in mols]))
    canon_smiles = []
    for s in smiles:
        try:
            canon_smiles.append(Chem.CanonSmiles(s))
        except:
            print(s)
    return canon_smiles

""" Check if two molecules in a set (combo: (m1,m2)) are equal to each other
    Checks through seeing if both have a substruct match with each other
    Input: set of mol objects (m1,m2)
    Output: True or False, depending on if they are equal (True) or not (False)
"""
def find_overlap(combo):
    m1, m2 = combo
    return m1.HasSubstructMatch(m2) and m2.HasSubstructMatch(m1)

def main():
    ### Number of fragments in each compound ###
    # #Read in fragments (Archaea & bacteria as a test)
    # fp = "Biology/Data/adenine_fragments.p"
    # #a_df = pd.read_csv(fp)
    # a_frags = pickle.load(open(fp, "rb")) #a_df["Frags"].tolist() #Smarts strings found in archaea
    # a_frags = [f for (m, f) in a_frags] #Separate mol files from fragments
    # a_mols = convert_to_mol_smarts(a_frags)
    # a_smiles = get_canonical_smiles(a_mols) #Convert into smiles (for later comparison)
    # print("Found", len(a_smiles), "adenine fragments")

    #Archaea smiles
    fp = "Biology/Data/archaea_cpds.csv"
    a_cpds = get_smiles(fp) #Smiles strings of archaea cpds
    a_cpds = [c for c in a_cpds if str(c) != 'nan'] #Remove nans
    a_cpds = convert_to_mol_smiles(a_cpds)
    print("Found", len(a_cpds), "archaea compounds")

    # #Find number of fragments in each compounds
    # frag_count = []
    # frags_in_all_cpds = []
    # for cpd in tqdm(a_cpds):
    #     count, frags_in_cpd = find_frags_in_cpds(cpd, a_mols)
    #     frag_count.append(count)
    #     frags_in_all_cpds.append(frags_in_cpd)
    #
    # #Basic descriptive statistics
    # print(stats.describe(frag_count))
    # pickle.dump(frag_count, open("Biology/Data/bacteria_fragment_counts.p", "wb"))
    # pickle.dump(frags_in_all_cpds, open("Biology/Data/bacteria_listOfAllFragsPerCpd.p", "wb"))
    #
    # ## Network of Compounds ###
    # frags_in_all_cpds = pickle.load(open("Biology/Data/bacteria_listOfAllFragsPerCpd.p", "rb"))
    #
    # #pickle.dump(create_graph(frags_in_all_cpds), open("Biology/Data/bacteria_frag_network.p", "wb"))
    # G = create_graph(frags_in_all_cpds)
    # nx.write_gexf(G, "Biology/Data/bacteria_frag_network.gexf")

    ### Molecular Assembly fragment comparison ###
    #Also check out https://stackoverflow.com/questions/51681659/how-to-use-rdkit-to-calculte-molecular-fingerprint-and-similarity-of-a-list-of-s for Tanimoto similarity values
    df = pd.DataFrame(columns=["Name","Total Fragments","Overlapping Fragment Count","Overlapping Fragments"])

    #for fp in [fp for fp in os.listdir("Other/Assembly_Fragments/") if ".txt" in fp]:

    fp = "Other/Assembly_Fragments/1Adenine_histWhole.txt" #adenine test
    assembly_mols = get_assembly_frags(fp)

    all_combinations = itertools.product(a_mols, assembly_mols)

    # print(assembly_smiles)
    #Find overlap - using dual substruct matching
    overlap_count = 0
    overlapping_frags = []
    for combo in all_combinations:
        if find_overlap(combo):
            overlap_count += 1
            m1, m2 = combo
            overlapping_frags.append(Chem.MolToSmiles(m2))

    #Find overlap between the two
    print(fp)
    print("Number of assembly fragments:", len(assembly_mols))
    print("Number of archaea fragments:", len(a_mols))
    print("Number of overlapping fragments:", overlap_count)
    print("Overlapping fragments: ", overlapping_frags)


    #new_entry = {"Name":fp, "Total Fragments:": len(assembly_smiles), "Overlapping Fragment Count":len(overlap), "Overlapping Fragments":overlap}
    df.loc[len(df)] = ["Adenine",len(assembly_mols),overlap_count,overlapping_frags]

    print(df.head())
    df.to_csv("Other/Assembly_Fragments/adenine_overlap.csv")

if __name__ == "__main__":
    main()
