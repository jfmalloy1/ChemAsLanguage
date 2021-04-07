from multiprocessing import Pool
from rdkit import Chem
from rdkit.Chem import MCS
import pickle
import time
from tqdm import tqdm
from itertools import combinations
from random import sample


""" Maximal common substring algorithm
    Input: mol objects of two compounds
    Output: mol object of largest common substring
"""
def fmcstimeout(p,q):
    return MCS.FindMCS([p,q], timeout=0.1).smarts

""" Wrapper for MCS function
    Input: set of two mols (c)
    Output: fragment of the two mol objects
"""
def findmcs(c):
    p,q = c
    #@timeout(2)
    try:
        return fmcstimeout(p,q)
    except:
        print("MCS of", p, "and", q, "timed out.")
        pass

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

def main():
    # ### AGAVE TEST ###
    # agave_test()

    # #Test: Time Fragments for Earth atmosphere SMILES strings
    # cpd_smiles = open("Other/Earth_atmosphere_SMILES.txt", "rb").readlines()

    ## KEGG smiles (~18k smiles strings, therefore 158 million combinations)
    kegg_smiles = open("Biology/Data/kegg_smiles.txt")
    cpd_smiles = kegg_smiles.readlines()

    cpd_mols = [Chem.MolFromSmiles(smi.strip()) for smi in cpd_smiles]
    cpd_mols = [m for m in cpd_mols if m != None]
    print("Retieved",len(cpd_mols),"random molecules in sample")

    ### BACKUP METHOD - in case paralleization doesn't work in time ###
    #random_fragment_generation(cpd_mols, 10)

    # # start = time.time()
    # # frags = []
    # # for s in (fragments(cpd_mols)): #| loadSmarts(sys.argv[2])):
    # #     try:
    # #         frags.append((Chem.MolFromSmarts(s),s))
    # #     except:
    # #         print("AAAGGGHHH",s,"failed")
    # # frags=UniqSmarts(frags)
    # # print("Found", len(frags), "many fragments")
    # # print("Time:", time.time() - start)
    #
    #

    ### ADENINE TEST (FOR ERNEST) ###
    #adenine_fragments("C1=NC2=NC=NC(=C2N1)N", cpd_mols)

    ### PARALLEL FRAGMENT GENERTION ###
    pool = Pool(processes=8)
    ## Sample from kegg smiles - from 1k to 5k (initially)
    for i in range(10):
        print("Iteration", i)
        for j in [1000]:#, 2000, 3000, 4000, 5000]:
            start = time.time()
            mols = sample(cpd_mols, i) #Sample i compounds to make fragments
            print(len(mols), "compounds being analyzed")
            cpd_combinations = combinations(mols, 2)
            frag_smarts = pool.map(findmcs, cpd_combinations)
            frags = []
            for s in (frag_smarts): #| loadSmarts(sys.argv[2])):
                try:
                    frags.append((Chem.MolFromSmarts(s),s))
                except:
                    pass
            frags=UniqSmarts(frags)
            print("Found", len(frags), "many fragments in", i, "compounds")
            print("Time:", time.time() - start)
            #
            pickle.dump(frags, open("Biology/Data/KEGG_fragments_" + str(i) + "cpds_iter" + str(i) + ".p", "wb"))

if __name__ == "__main__":
    main()
