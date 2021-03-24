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
    fp = "../br08001.json"
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

""" Find & label all unique domain compounds
    Input: none (assumes archaea/bacteria/eukaraya_cpds.csv as files in working directory)
    Output: Dictionary of unique domain cpds, labeled with specific domains (archaea/bacteria/eukarya)
"""
def load_domain_cpds():
    #Dataframe of each set of cpds
    a_df = pd.read_csv("archaea_cpds.csv")
    b_df = pd.read_csv("bacteria_cpds.csv")
    e_df = pd.read_csv("eukarya_cpds.csv")

    #Uniqe cpds of each
    a_nobacteria = a_df[~a_df.compounds.isin(b_df.compounds)]
    a_unique = a_nobacteria[~a_nobacteria.compounds.isin(e_df.compounds)]
    print(len(a_unique))

    b_noarchaea = b_df[~b_df.compounds.isin(a_df.compounds)]
    b_unique = b_noarchaea[~b_noarchaea.compounds.isin(e_df.compounds)]
    print(len(b_unique))

    e_nobacteria = e_df[~e_df.compounds.isin(b_df.compounds)]
    e_unique = e_nobacteria[~e_nobacteria.compounds.isin(a_df.compounds)]
    print(len(e_unique))

    #Make a dictionary with each compound, labeled accordingly
    domain_cpds = {}
    for index, row in a_unique.iterrows():
        domain_cpds[row["compounds"][:6] + ".mol"] = "Archaea"
    for index, row in b_unique.iterrows():
        domain_cpds[row["compounds"][:6] + ".mol"] = "Bacteria"
    for index, row in e_unique.iterrows():
        domain_cpds[row["compounds"][:6] + ".mol"] = "Eukarya"

    return domain_cpds


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

""" Find & add all fragment vectors within a compound
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
            try:
                vs.append(word2vec[f])
            except:
                pass
        #Only sum vectors if vectors were present in the compound
        if vs:
            vectors.append(np.sum(vs, axis=0).astype("float64"))

    return vectors

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

def main():
    #Load w2v model, kegg dataframe, and all fragments
    word2vec = load_w2v()
    kegg_df = pd.read_csv("../kegg_data.csv")
    frags = get_chem_fragments("../frags.txt")

    #Find unique domain compounds, load them into a dictionary
    domain_cpds = load_domain_cpds()

    # #Compound class dictionary
    # cpd_classes = find_cpd_classes()

    #Find vectors associated with all domain cpds
    full_df = pd.DataFrame(columns=np.arange(0,100))
    for d in ["Archaea", "Bacteria", "Eukarya"]:
        full_df = full_df.append(get_class_dataframe(kegg_df, word2vec, frags, d, domain_cpds))

    print(full_df)

    # #Add compound classes
    # for class_label in list(set(cpd_classes.values())):
    #     full_df = full_df.append(get_class_dataframe(kegg_df, word2vec, frags, class_label, cpd_classes))
    #
    ## RANDOM CPDS ##
    #Find 1000 of them, to distinguish between 438 LUCA cpds
    # rand_cpds = find_random_SMILES(kegg_df, 100, list(luca_cpds.keys()))
    # rand_frags = find_frags_within_SMILES(rand_cpds, frags)
    # rand_vectors = add_frag_vectors(word2vec, rand_frags)
    # rand_df = pd.DataFrame(rand_vectors)
    # rand_df["label"] = ["Random"] * len(rand_df)
    # #print("Number of random vectors:", len(rand_df))

    # #Combine luca and random dataframes
    # full_df = full_df.append(rand_df)

    #Run TSNE - 2 classes (LUCA and random)
    TSNE_visual(full_df, len(full_df["label"].unique()))


if __name__ == "__main__":
    main()
