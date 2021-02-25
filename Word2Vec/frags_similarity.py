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

""" Find & label all 438 LUCA compounds
    Input: none (assumes Luca_ec_compound_unique_list_clean.json as a file in working directory)
    Output: Dictionary of LUCA cpds labeled with "LUCA"
"""
def load_luca_cpds():
    fp = "LUCA/Luca_ec_compound_unique_list_clean.json"
    with open(fp) as json_file:
        luca_json = json.load(json_file)

    #Label all luca compounds with "LUCA"
    luca_cpds = {}
    for cpd in luca_json["Luca"]:
        luca_cpds[cpd[:6] + ".mol"] = "LUCA"

    return luca_cpds

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

""" Find the vector associated with a fragment
    Goal is to have a single vector for each fragment
    Inputs: trained word2vec model, fragment to find
    Outputs: one vector per fragment
"""
def find_frag_vector(word2vec, frag):
    return word2vec[frag]

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

    plt.figure(figsize=(16,10))
    sns.scatterplot(
        x="tsne-2d-one", y="tsne-2d-two",
        hue="label",
        palette=sns.color_palette("hls", n_categories),
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

    #Find and label all fragments associated with a cpd class
    class_vectors = []
    for cpd in class_frags:
        for f in cpd:
            class_vectors.append(find_frag_vector(word2vec, f))

    class_df = pd.DataFrame(class_vectors)
    class_df["label"] = [class_label] * len(class_df)
    print("Number of", class_label, "compounds:", len(class_df))
    return class_df

def main():
    word2vec = load_w2v()
    kegg_df = pd.read_csv("kegg_data.csv")

    #Load LUCA compounds into a dictionary
    luca_cpds = load_luca_cpds()

    #Get chemical fragment list
    frags = get_chem_fragments("frags.txt")

    #Compound class dictionary
    cpd_classes = find_cpd_classes()

    #Create empty, 100 column dataframe to hold vectors
    full_df = pd.DataFrame(columns=np.arange(0,100))

    # for class_label in list(set(cpd_classes.values())):
    #     full_df = full_df.append(get_class_dataframe(kegg_df, word2vec, frags, class_label, cpd_classes))

    ## RANDOM CPDS ##
    rand_cpds = find_random_SMILES(kegg_df, 100, list(cpd_classes.keys()))
    rand_frags = find_frags_within_SMILES(rand_cpds, frags)
    rand_vectors = []
    for cpd in rand_frags:
        for f in cpd:
            rand_vectors.append(find_frag_vector(word2vec, f))
    rand_df = pd.DataFrame(rand_vectors)
    rand_df["label"] = ["Random"] * len(rand_df)
    print("Number of random fragments:", len(rand_df))

    #Combine aa and random dataframes
    full_df = full_df.append(rand_df)

    #Add LUCA vectors
    full_df = full_df.append(get_class_dataframe(kegg_df, word2vec, frags, "LUCA", luca_cpds))

    #Run TSNE - number of classes is the number of cpd classes + 1 (random compounds)
    #TSNE_visual(full_df, len(list(set(cpd_classes.values()))) + 1) #Color palette involves classes
    TSNE_visual(full_df, full_df["label"].nunique())


if __name__ == "__main__":
    main()
