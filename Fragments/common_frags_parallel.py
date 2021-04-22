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


""" Tests if two fragments are the same """
def sameMolecule(a,b):
    def same_or_timeout(a,b):
        if a[1] == b[1] : return True #string compare
        if a[0].GetNumAtoms() != b[0].GetNumAtoms() : return False
        if a[0].GetNumBonds() != b[0].GetNumBonds() : return False
        return a[0].HasSubstructMatch(b[0]) and b[0].HasSubstructMatch(a[0])

    try :
        res = same_or_timeout(a,b)
    except:
        res = False

    return res

""" Ensures unique fragments
    Input: set of fragments
    Output: unique set of fragments
"""
def UniqSmarts(frags):
    result = set()
    frags = set((g,t) for (g,t) in frags if g != None) # was successfully converted to smarts
    while frags:
        f,s = frags.pop()
        result.add((f,s))
        frags = set((g,t) for (g,t) in frags
                          if not (sameMolecule((f,s),(g,t)))) # check if they are the same or not. check also on the string to speed up
    return result

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
def find_unique_frags(input_fp, output_fp):
    start = time.time()
    print("Analyzing:", input_fp)
    frag_smarts = pickle.load(open(input_fp, "rb"))
    print("Time to load:", time.time() - start)
    #Remove fragments with the same smarts strings
    print("Found", len(frag_smarts), "original fragments")
    frag_smarts = list(set(frag_smarts))
    print("Found", len(frag_smarts), "fragments after removing duplicate SMARTS")
    frags = []
    for s in tqdm((frag_smarts)): #| loadSmarts(sys.argv[2])):
        try:
            frags.append((Chem.MolFromSmarts(s),s))
        except:
            pass
    frags=UniqSmarts(frags)
    print("Found", len(frags), "many unique fragments")
    print("Time:", time.time() - start)

    pickle.dump(frags, open(output_fp, "wb"))

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

    kegg_size = len(kegg_mols)
    del kegg_mols
    for i in range(10):
        reaxys_mols = sample_Reaxys(df, kegg_size)

        ### ADENINE TEST (FOR ERNEST) ###
        #adenine_fragments("C1=NC2=NC=NC(=C2N1)N", cpd_mols)

        ### PARALLEL FRAGMENT GENERTION ###
        pool = Pool(processes=16)

        generate_fragments(pool, reaxys_mols, "Technology/Data/Reaxys_fragments_keggSize_" + str(i) + ".p")

        ### FIND UNIQUE FRAGMENTS ###
        find_unique_frags("Technology/Data/Reaxys_fragments_keggSize_" + str(i) + ".p", "Technology/Data/Reaxys_fragments_keggSize_0" + str(i) + "_unique.p")

if __name__ == "__main__":
    main()
