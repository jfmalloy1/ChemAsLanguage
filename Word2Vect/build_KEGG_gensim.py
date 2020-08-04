""" Use for gensim W2V model
    Goal - output a list of fragments/functional groups, in order of appeareance within molecules
"""
import pandas as pd
from rdkit import Chem
import numpy as np
import re
import itertools
import random
import os
from gensim.models import Word2Vec
from gensim.test.utils import get_tmpfile
from gensim.models import KeyedVectors

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

""" Return the mol representation of a given Smiles string
    Input: Smiles string of a compound
    Output: mol representation
"""
def get_mol_representation(smiles):
    return Chem.MolFromSmiles(smiles)

""" Return an ordered list of all fragments present in the given compounds
    Input: Dataframe containing smiles strings of compounds
    Output: ordered list of chemical fragments
"""
def get_fragment_list(df, fragments):
    #Master list of words
    words = []
    #Remove old line-sentence format
    os.remove("frags_linesentence.txt")

    #Loop through all molecules
    for index, row in df.iterrows():
        print(index)
        #store fragments present in each molecule
        mol_words = []
        #Turn smiles into mol format
        mol = get_mol_representation(row["Original SMILES"])
        #find fragments present
        for f in fragments:
            try:
                if mol.HasSubstructMatch(Chem.MolFromSmarts(f)):
                    mol_words.append(f)
            except:
                continue
        #Randomize order of fragments (because what is order in a compound?)
        random.shuffle(mol_words)
        words.append(mol_words)

        #Make a file of the fragements within each molecule - easier to train gensim model on
        for w in mol_words:
            print(w, end=" ", file=open("frags_linesentence.txt", "a"))
        print(file=open("frags_linesentence.txt", "a"))

def main():
    #Read in KEGG data
    df = get_KEGG_smiles("kegg_data.csv")
    #Read in chem fragments
    fragments = get_chem_fragments("frags.txt")

    ### TEST on a sample of KEGG ###
    ordered_frags = get_fragment_list(df, fragments)
    print(ordered_frags)
    # for frag in ordered_frags:
    #     print(frag)

    # #Build Word2Vec model on ordered fragments
    word2vec = Word2Vec(corpus_file = "frags_linesentence.txt", min_count=2)

    ## MODEL TESTING ##
    v1 = word2vec.wv["[#6]-[#6]"]
    print(v1)

    sim_frags = word2vec.wv.most_similar("[#6]-[#6]")
    print(sim_frags)

    #Save trained model
    word_vectors = word2vec.wv
    #fname = get_tmpfile("vectors_01.kv")
    word_vectors.save("vectors_fullKEGG.kv")


if __name__ == "__main__":
    main()
