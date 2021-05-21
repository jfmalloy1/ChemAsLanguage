from multiprocessing import Pool
from rdkit import Chem
from rdkit.Chem import MCS
from rdkit import RDLogger
import pandas as pd
import pickle
import time
from tqdm import tqdm
from itertools import combinations
from itertools import repeat
from random import sample
import os


""" Maximal common substring algorithm
    Input: mol objects of two compounds
    Output: mol object of largest common substring
"""
def fmcstimeout(p,q):
    return MCS.FindMCS([p,q], timeout=0.01).smarts

""" Wrapper for MCS function
    Input: set of two mols (c)
    Output: fragment of the two mol objects
"""
def findmcs(c):
    p,q = c
    #@timeout(2)
    # try:
    return fmcstimeout(p,q)
    # except:
    #     print("MCS of", p, "and", q, "timed out.")
    #     pass

""" Finds the fragments within a set of mol objects
    Input: A set of rdkit mol objects
    Output: set of fragments (in mol object) associated with the OG set of mols
    NOTE: CHANGED FROM p,q in combinations() TO c in combinations() TO ACCOUNT FOR PARALLELIZATION
"""
def fragments(mols):
    s=[]
    for c in combinations(mols,2):
        try:
            s.append(findmcs(c))
        except:
            pass

    return set(s)

""" Tests if two mol objects (generated from fragments) are the same """
def same_or_timeout(c):
    a,b = c
    if a.HasSubstructMatch(b) and b.HasSubstructMatch(a):
        return (Chem.MolToSmarts(a),Chem.MolToSmarts(b))
    else:
        return False

""" RDKit test on ASU agave cluster
    Input: None
    Output: None
"""
def agave_test():
    m = Chem.MolFromSmiles("Cc1ccccc1")
    print(Chem.MolToSmiles(m))

""" Assocaited with Ernest's Adenine test - makes specific combinations with adenine
    Input: Adenine smiles string, all KEGG cpd mol objects
    Output: Set of combinations of [(Adenine, C1), (Adenine, C2)], etc...
"""
def make_combinations(adenine_smiles, cpd_mols):
    adenine_mol = Chem.MolFromSmiles(adenine_smiles)
    combos = []
    for m in cpd_mols:
        combos.append(set([adenine_mol, m]))
    return combos

""" Test for Ernest - find all fragments associated with Adenine
    Input: Adenine smiles string, all KEGG compound mol objects
    Output: Set of fragments associated with Adenine
"""
def adenine_fragments(adenine_smiles, cpd_mols):
    start = time.time()
    pool = Pool(processes=4)

    adenine_combinations = make_combinations(adenine_smiles, cpd_mols)
    frag_smarts = pool.map(findmcs, adenine_combinations)
    frags = []
    for s in (frag_smarts): #| loadSmarts(sys.argv[2])):
        try:
            frags.append((Chem.MolFromSmarts(s),s))
        except:
            pass
    frags=UniqSmarts(frags)
    print(frags)
    print("Found", len(frags), "many fragments")
    print("Time:", time.time() - start)
    #
    pickle.dump(frags, open("Biology/Data/adenine_fragments.p", "wb"))

""" BACKUP METHOD: compares a specified number of pairwise molecules to generate Fragments
    Input: list of mol objects, number of pairwise comparisons (n)
    Output: fragments saved to a pickle file
"""
def random_fragment_generation(cpd_mols, n):
    for i in range(n):
        mols = sample(cpd_mols, 2)
        frag = fmcstimeout(mols[0], mols[1])
        print(frag)

""" Read in KEGG mol objects
    Input: None (assumes list of all KEGG smiles strings is in "Biology/Data/kegg_smiles.txt")
    Output: List of RDKit mol objects of all of KEGG
"""
def read_KEGG_mols():
    ## KEGG smiles (~18k smiles strings, therefore 158 million combinations)
    kegg_smiles = open("Biology/Data/kegg_smiles.txt")
    cpd_smiles = kegg_smiles.readlines()

    cpd_mols = [Chem.MolFromSmiles(smi.strip()) for smi in cpd_smiles]
    cpd_mols = [m for m in cpd_mols if m != None]
    print("Retieved",len(cpd_mols),"random molecules in sample")
    return cpd_mols

""" Generate fragments in parallel
    Input: pool (parallelization), mols (list of mol objects to be used), output_fp (filepath to place fragments)
    Output: fragments in a pickled file
"""
def generate_fragments(pool, mols, output_fp):
    start = time.time()
    print(len(mols), "compounds being analyzed")
    cpd_combinations = combinations(mols, 2) #Generate combinations
    frag_smarts = pool.map(findmcs, cpd_combinations) #Find fragments in parallel
    print("Found", len(frag_smarts), "possible fragments in", time.time() - start, "seconds")
    #Save smarts fragments
    pickle.dump(frag_smarts, open(output_fp, "wb"))

""" Find the unique fragments in a set
    Input: input filepath (assumes this is a pickled list of fragments in SMARTS form), output filepath
    Output: List of unique fragments in pickled form
"""
def find_unique_frags(pool, input_fp, output_fp):
    start = time.time()
    print("Analyzing:", input_fp)
    frag_smarts = pickle.load(open(input_fp, "rb"))
    print("Time to load:", time.time() - start)

    #Remove fragments with the same smarts strings
    print("Found", len(frag_smarts), "original fragments")
    frag_smarts = list(set(frag_smarts))
    print("Found", len(frag_smarts), "fragments after removing duplicate SMARTS")

    frag_mols = []
    for s in frag_smarts: #| loadSmarts(sys.argv[2])):
        try:
            frag_mols.append(Chem.MolFromSmarts(s))
        except:
            pass

    cpd_combinations = combinations(frag_mols, 2) #Generate combinations

    #Parallel unique function - chunksize=1000 based on 10 ms best practice (1 iter=0.01ms)
    non_unique_frags = list(pool.imap(same_or_timeout, cpd_combinations, chunksize = 1000))
    non_unique_frags = list(filter(bool,non_unique_frags))

    print("Found", len(non_unique_frags), "nonunique fragments")

    #For each nonunique set, remove one of the nonunique fragments (leave the other)
    for c in non_unique_frags:
        s1, s2 = c
        try:
            frag_smarts.remove(s1)
        except:
            try:
                frag_smarts.remove(s2)
            except:
                #If neither were in the fragment list, problem solved
                pass

    print("Found", len(frag_smarts), "many unique fragments")
    print("Time:", time.time() - start)

    pickle.dump(frag_smarts, open(output_fp, "wb"))

""" Creates a dataframe from the Reaxys subset number n between 1-10
    Input: the number of the specific section of the database
    Output: dataframe containing all of Reaxys smiles strings
"""
def read_cpds(n):
    start = time.time()
    df = pd.read_csv("~/Shared/Reaxys/cpd_"+n+".csv", header=None, usecols=[31])

    df.columns=["smiles"]
    end = time.time()
    print("Time for " + n + ":", end-start)
    return df

""" Read and sample Reaxys, given a specific sample size. Assumes Reaxys is in a specific directory.
    Input: Reaxys dataframe, sample size (s)
    Output: list of RDKit mol objects
"""
def sample_Reaxys(df, s):
    #Remove rdkit warnings
    RDLogger.DisableLog('rdApp.*')

    #Sample given sample size s
    smiles = df.sample(frac=s/len(df.index))["smiles"].tolist()
    smiles = list(map(str, smiles))

    #Convert all sampled smiles strings into mols
    mols = [Chem.MolFromSmiles(smi.strip()) for smi in smiles]
    mols = [m for m in mols if m != None]
    print("Retieved",len(mols),"random molecules")
    return mols

def main():
    # ### AGAVE TEST ###
    # agave_test()

    # #Test: Time Fragments for Earth atmosphere SMILES strings
    # cpd_smiles = open("Other/Earth_atmosphere_SMILES.txt", "rb").readlines()

    ### Get KEGG Mol Objects ###
    kegg_mols = read_KEGG_mols()

    ### Get Reaxys Mol Objects ###
    #Read in full reaction database
    df = pd.DataFrame()
    for i in range(1,11):
        df = df.append(read_cpds(str(i)), ignore_index=True)
        print("Done with subset", i, "...")
        print("Df size", len(df.index))

    ### ADENINE TEST (FOR ERNEST) ###
    #adenine_fragments("C1=NC2=NC=NC(=C2N1)N", cpd_mols)

    ### PARALLEL FRAGMENT GENERTION ###
    pool = Pool(processes=16)
    RDLogger.DisableLog('rdApp.*')

    #kegg_size = len(kegg_mols)
    for i in range(10):
        print("Analyzing sample", i)
        fp = "Technology/Data/Reaxys_1000_Samples/"
        reaxys_mols = sample_Reaxys(df, 1000)

        #Save mols for future occurrence testing
        pickle.dump(reaxys_mols, open(fp + "sample_" + str(i) + "_ReaxysMols.p", "wb"))

        generate_fragments(pool, reaxys_mols, fp + "sample_" + str(i) + "frags.p")

        ### FIND UNIQUE FRAGMENTS ###
        find_unique_frags(pool, fp + "sample_" + str(i) + "frags.p", fp + "sample_" + str(i) + "frags_unique.p")
        print()

if __name__ == "__main__":
    main()
