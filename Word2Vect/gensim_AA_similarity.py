import pandas as pd
from rdkit import Chem
import numpy as np
from gensim.models import Word2Vec
from gensim.test.utils import get_tmpfile
from gensim.models import KeyedVectors

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
def find_frags_within_SMILES(aa_smiles, frags):
    aa_frags = []
    for aa in aa_smiles:
        #Turn AA into mol file
        mol = Chem.MolFromSmiles(aa)

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
        aa_frags.append(list(set(individual_frags)))

    return aa_frags

""" Find all SMILES sequences for amino acids
    Input: None - list of AA ids are from KEGG
    Output: list of all SMILES strings for AAs
"""
def find_AA_SMILES():
    ## AA strings ##
    AA_ids = ["C00037", "C00041", "C00183", "C00123", "C00407", "C00049", "C00152",
        "C00025", "C00064", "C00065", "C00188", "C00073", "C00097", "C00047", "C00062",
        "C00135", "C00148", "C00079", "C00082", "C00078"]
    #Add .mol onto every AA id
    for i in range(len(AA_ids)):
        AA_ids[i] = AA_ids[i] + ".mol"

    kegg_df = pd.read_csv("kegg_data.csv")
    #Find all AAs
    aa_df = kegg_df[kegg_df["MOL file"].isin(AA_ids)]
    return aa_df["Original SMILES"].tolist()

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
        vectors.append(np.sum(vs, axis=0))

    return vectors

def main():
    word2vec = load_w2v()

    # ## MODEL TEST ##
    # #Find vector & similarities of C-C bond
    # v1 = word2vec.wv["[#6]-[#6]"]
    # print(v1)
    # sim_frags = word2vec.wv.most_similar("[#6]-[#6]")
    # print(sim_frags)

    #Get chemical fragment list
    frags = get_chem_fragments("frags.txt")

    #Get amino acid smiles strings
    aas = find_AA_SMILES()

    #Find all fragments within each amino acid
    aa_frags = find_frags_within_SMILES(aas, frags)

    #Find vector sum of each AA
    aa_vectors = add_frag_vectors(word2vec, aa_frags)
    for aa_v in aa_vectors:
        print(aa_v)
        print("-----------------------------")


if __name__ == "__main__":
    main()
