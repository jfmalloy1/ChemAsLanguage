import pandas as pd
from rdkit import Chem
import numpy as np
import json
from gensim.models import Word2Vec
from gensim.test.utils import get_tmpfile
from gensim.models import KeyedVectors
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import re

""" Load trained Word2Vec model - Gensim on full KEGG
    Input: None
    Output: trained word2vec model
"""
def load_w2v():
    return KeyedVectors.load("../vectors_fullKEGG.kv", mmap="r")

""" Return a list of all chemical fragments in SMARTS form """
""" Input: filepath (fp) to txt file """
""" Output: list of SMARTS fragments """
def get_chem_fragments(fp):
    fragments = []
    with open(fp, "r") as f:
        for line in f:
            fragments.append(line.strip())

    return fragments

""" Find the fragments within a list of smiles strings
    Input: SMILES string
    Output: List of lists of fragments within smiles strings
"""
def find_frags_within_SMILES(cpd_list, smiles, frags):
    cpd_frags = []
    removed_cpds = []
    i = 0
    for smi in smiles:
        #Turn AA into mol file
        mol = Chem.MolFromSmiles(smi)

        #Loop through all fragments to find occurances within AAs
        individual_frags = []
        for f in frags:
            try:
                #If a fragment is found in an AA, add it to the individual frags list
                if mol.HasSubstructMatch(Chem.MolFromSmarts(f)):
                    individual_frags.append(f)
            except:
                removed_cpds.append(cpd_list[i])
                pass
        #Add each individual AA to AA frags - remove
        cpd_frags.append(list(set(individual_frags)))
        i += 1

    return cpd_frags, [x for x in cpd_list if x not in removed_cpds]

""" Find all SMILES sequences for a random subset of KEGG
    Input: kegg dataframe, number of samples to be collected
    Output: a list of smiles strings from the random sample
"""
def find_random_SMILES(kegg_df, n_samples, cpds_to_ignore):
    #Remove all compounds which are being classified
    kegg_df = kegg_df[~kegg_df["MOL file"].isin(cpds_to_ignore)]
    #Remove all empty SMILES values
    kegg_df = kegg_df.dropna(subset=["Original SMILES"])
    kegg_df = kegg_df[kegg_df["Original SMILES"] != ""]
    #Randomly sample KEGG
    sub_df = kegg_df.sample(n_samples)

    #Return smiles strings
    return sub_df["Original SMILES"].tolist(), sub_df["MOL file"].tolist()

""" Find & add all fragment vectors within a compound
    Goal is to have a single vector for each compound
    Inputs: trained word2vec model, list of lists of fragments within amino acids
    Outputs: one vector (sum of all fragment vectors) per amino acid
"""
def add_frag_vectors(cpd_list, word2vec, frags):
    vectors = []
    removed_cpds = []
    i = 0
    #loop through amino acids
    for cpd in frags:
        vs = []
        #Loop through fragments within each amino acid, add vectors to a list
        for f in cpd:
            try:
                vs.append(word2vec[f])
            except:
                pass
        #Only sum vectors if vectors were present in the compound
        if vs:
            vectors.append(np.sum(vs, axis=0).astype("float64"))
        else:
            removed_cpds.append(cpd_list[i])

        #Ensure the correct compound gets removed
        i+=1

    return vectors, [x for x in cpd_list if x not in removed_cpds]

""" Run TSNE visualization
    Input: dataframe of compoud vectors (df["label"] is the compound label)
    Output: Visualization of the trained vectors
"""
def TSNE_visual(df, n_categories):
    #find values to pass to TSNE
    data_values = df[list(range(0,100))].values

    tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
    tsne_results = tsne.fit_transform(data_values)

    df["tsne-2d-one"] = tsne_results[:,0]
    df["tsne-2d-two"] = tsne_results[:,1]

    pal = sns.color_palette("hls", n_categories)

    plt.figure(figsize=(16,10))
    sns.scatterplot(
        x="tsne-2d-one", y="tsne-2d-two",
        hue="label",
        palette=sns.color_palette(palette=pal),
        data=df,
        legend="full"
    )
    plt.show()

""" Find all compounds associated with a particular class of compounds within KEGG
    Input: dataframe of KEGG data, trained W2V model, fragment list, label of class to search for
    Output: dataframe of vectors associated with a particular class
"""
def get_class_dataframe(kegg_df, word2vec, frags, class_label, cpd_classes):
    #Find all compound IDs associated with a particular label
    cpd_ids = [k for k,v in cpd_classes.items() if v == class_label]
    cpd_smiles = kegg_df[kegg_df["MOL file"].isin(cpd_ids)]["Original SMILES"].tolist()

    class_frags = find_frags_within_SMILES(cpd_smiles, frags)
    vectors = add_frag_vectors(word2vec, class_frags)
    class_df = pd.DataFrame(vectors)
    class_df["label"] = [class_label] * len(class_df)
    print("Number of", class_label, "compounds:", len(class_df))
    return class_df

""" Builds a KEGG network, finds the central nodes, calculates distances between all nodes
    Input: None (assumes newKEGG_reactionEdges.json exists within the current directory
    Output: networkx graph of KEGG (unipartite, consisting only of compounds), the central node of the network, distances between all cpds
"""
def KEGG_network():
    #Load KEGG json file
    f = open("newKEGG_reactionEdges.json", "r")
    kegg = json.load(f)

    #Find all reaction-compound pairs
    rxn_list = []
    cpd_list = []
    cpd_rxn_pairs = []

    #loop through both products and substrates
    for option in ["products", "substrates"]:
      for rxn in kegg[option]:
        #add each reaction to a master list
        rxn_list.append(rxn)
        for cpd in kegg["products"][rxn]:
          #add each compound to a master list
          cpd_list.append(cpd)
          #create a tuple of each cpd_rxn pair, add them to a master list
          cpd_rxn_pairs.append(tuple([cpd, rxn]))

    #remove duplicates of reactions and compounds
    rxn_list = list(set(rxn_list))
    cpd_list = list(set(cpd_list))

    #Create a bipartite graph using reactions, compounds, and the cpd/rxn pair
    KEGG_graph = nx.Graph()
    KEGG_graph.add_nodes_from(rxn_list, bipartite=0)
    KEGG_graph.add_nodes_from(cpd_list, bipartite=1)
    KEGG_graph.add_edges_from(cpd_rxn_pairs)

    #Create a project of only compounds
    KEGG_cpd_graph = nx.bipartite.projected_graph(KEGG_graph, cpd_list)

    #Find the central node of the largest connected component
    lcc = max(nx.connected_components(KEGG_cpd_graph), key=len)
    lcc_graph = KEGG_cpd_graph.subgraph(lcc)
    ## CENTER(s) ##
    centers = ['C00006', 'C00014', 'C00025', 'C00001', 'C00011']

    #Calculate distances between all nodes
    distances = dict(nx.all_pairs_shortest_path_length(KEGG_cpd_graph))

    return KEGG_cpd_graph, centers, distances

""" Find the maximum distance between a given compound and the centers of the graph
    Input: Compound, centers of the largest connected component within the graph
    Output: Distance (int), or "NC" if not connected
"""
def find_distance(cpd, centers, distances):
    d = []
    #Find the distance between the "centers" of the largest connected component
    for c in centers:
        try:
            d.append(distances[c][cpd])
        except:
            pass

    #If the random compound is not in the largest connected component, labeld "NC" (not connected)
    if not d:
        return "NC"
    #Otherwise, label with the max distance from the center
    else:
        return str(max(d))

def main():
    #Load w2v model, kegg dataframe, and all fragments
    word2vec = load_w2v()
    kegg_df = pd.read_csv("../kegg_data.csv")
    frags = get_chem_fragments("../frags.txt")

    KEGG_cpd_graph, centers, distances = KEGG_network()

    ## RANDOM CPDS ##
    #Find 10 of them, ignoring no compounds (initially)
    rand_cpds, cpd_list = find_random_SMILES(kegg_df, 100, [])
    rand_frags, cpd_list = find_frags_within_SMILES(cpd_list, rand_cpds, frags)
    rand_vectors, cpd_list = add_frag_vectors(cpd_list, word2vec, rand_frags)
    rand_df = pd.DataFrame(rand_vectors)
    rand_df["Cpds"] = cpd_list

    #Label by max distance from central compound
    cpd_distance = []
    for index, row in rand_df.iterrows():
        cpd_distance.append(find_distance(re.sub(".mol", "", row["Cpds"]), centers, distances))

    rand_df["label"] = cpd_distance
    #print("Number of random vectors:", len(rand_df))

    #Run TSNE
    TSNE_visual(rand_df, len(rand_df["label"].unique()))


if __name__ == "__main__":
    main()
