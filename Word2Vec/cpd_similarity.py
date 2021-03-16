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
from sklearn import cluster

""" Load trained Word2Vec model - Gensim on full KEGG
    Input: None
    Output: trained word2vec model
"""
def load_w2v():
    return KeyedVectors.load("vectors_fullKEGG.kv", mmap="r")

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
def find_frags_within_SMILES(smiles, frags):
    cpd_frags = []
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
                pass
        #Add each individual AA to AA frags - remove
        cpd_frags.append(list(set(individual_frags)))

    return cpd_frags

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

""" Find all SMILES sequences for amino acids
    Input: None - list of AA ids are from KEGG
    Output: list of all SMILES strings for AAs
"""
def find_AA_SMILES(kegg_df):
    ## AA strings ##
    AA_ids = ["C00037", "C00041", "C00183", "C00123", "C00407", "C00049", "C00152",
        "C00025", "C00064", "C00065", "C00188", "C00073", "C00097", "C00047", "C00062",
        "C00135", "C00148", "C00079", "C00082", "C00078"]
    #Add .mol onto every AA id
    for i in range(len(AA_ids)):
        AA_ids[i] = AA_ids[i] + ".mol"

    #Find all AAs
    aa_df = kegg_df[kegg_df["MOL file"].isin(AA_ids)]
    return AA_ids, aa_df["Original SMILES"].tolist()

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
    return sub_df["Original SMILES"].tolist()

""" Find & add all fragment vectors within a single amino acid
    Goal is to have a single vector for each compound
    Inputs: trained word2vec model, list of lists of fragments within amino acids
    Outputs: one vector (sum of all fragment vectors) per amino acid
"""
def add_frag_vectors(word2vec, frags):
    vectors = []
    #loop through amino acids
    for cpd in frags:
        vs = []
        #Loop through fragments within each amino acid, add vectors to a list
        for f in cpd:
            vs.append(word2vec[f])
        #Only sum vectors if vectors were present in the compound
        if vs:
            vectors.append(np.sum(vs, axis=0).astype("float64"))

    return vectors

""" Run TSNE visualization
    Input: dataframe of compoud vectors (df["label"] is the compound label)
    Output: Visualization of the trained vectors
"""
def TSNE_visual(df, n_categories, classification):
    #find values to pass to TSNE
    data_values = df[list(range(0,100))].values

    tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
    tsne_results = tsne.fit_transform(data_values)

    df["tsne-2d-one"] = tsne_results[:,0]
    df["tsne-2d-two"] = tsne_results[:,1]

    plt.figure(figsize=(16,10))
    sns.scatterplot(
        x="tsne-2d-one", y="tsne-2d-two",
        hue=classification,
        palette=sns.color_palette("hls", n_categories),
        data=df
        #legend="full"
    )

    plt.xlabel("")
    plt.ylabel("")
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

""" K-means clustering of word2vec data
    Input: model (trained word2vec model), num_clusters (the number of clusters)
    Output: Output dictionary of compounds and cluster labels associated with them (to compare with W2V)
"""
def clustering(full_df, num_clusters):
    sub_df = full_df.drop(columns=["label"])
    X = sub_df.values.tolist()
    # print(model)
    # X = model[model.wv.vocab]
    # print(X)
    #set up clusters
    kmeans = cluster.KMeans(n_clusters = num_clusters)
    #fit Word2Vec data to cluster model
    kmeans.fit(X)

    #kmean_labels = kmeans.labels_
    full_df["kmean_labels"] = kmeans.labels_

    #Test accuratness of k-means (//TODO: talk to Yanbo about better ways of doing this?)
    for index, row in full_df.iterrows():
        print(row["label"], row["kmean_labels"])
    return full_df

def main():
    word2vec = load_w2v()
    kegg_df = pd.read_csv("kegg_data.csv")

    # ## MODEL TEST ##
    # #Find vector & similarities of C-C bond
    # v1 = word2vec.wv["[#6]-[#6]"]
    # print(v1)
    # sim_frags = word2vec.wv.most_similar("[#6]-[#6]")
    # print(sim_frags)

    #Get chemical fragment list
    frags = get_chem_fragments("frags.txt")

    #Compound class dictionary
    cpd_classes = find_cpd_classes()

    # ## AMINO ACIDS ##
    # #Get amino acid smiles strings
    # aa_ids, aas = find_AA_SMILES(kegg_df)
    # #Find all fragments within each amino acid
    # aa_frags = find_frags_within_SMILES(aas, frags)
    # #Find vector sum of each AA
    # aa_vectors = add_frag_vectors(word2vec, aa_frags)
    # #Create dataframe with AAs & vectors
    # aa_df = pd.DataFrame(aa_vectors)
    # aa_df["label"] = ["AA"] * len(aa_ids)
    # print("Number of AAs:", len(aa_df))

    full_df = pd.DataFrame(columns=np.arange(0,100))
    for class_label in list(set(cpd_classes.values())):
        full_df = full_df.append(get_class_dataframe(kegg_df, word2vec, frags, class_label, cpd_classes))

    ## RANDOM CPDS ##
    # rand_cpds = find_random_SMILES(kegg_df, 100, list(cpd_classes.keys()))
    # rand_frags = find_frags_within_SMILES(rand_cpds, frags)
    # rand_vectors = add_frag_vectors(word2vec, rand_frags)
    # rand_df = pd.DataFrame(rand_vectors)
    # rand_df["label"] = ["Random"] * len(rand_df)
    # print("Number of random compounds:", len(rand_df))

    # #Combine aa and random dataframes
    # full_df = full_df.append(rand_df)

    #Run TSNE - number of classes is the number of cpd classes + 1 (random compounds)
    TSNE_visual(full_df, len(list(set(cpd_classes.values()))), "label")# + 1)

    #K-means clustering
    #full_df = clustering(full_df, len(list(set(cpd_classes.values()))))
    #TSNE_visual(full_df, len(list(set(cpd_classes.values()))), "kmean_labels")


if __name__ == "__main__":
    main()
