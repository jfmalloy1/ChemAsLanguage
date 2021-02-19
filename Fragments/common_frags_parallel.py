from multiprocessing import Pool
from rdkit import Chem
from rdkit.Chem import MCS
import pickle
import time
from tqdm import tqdm
from itertools import combinations


""" Maximal common substring algorithm
    Input: mol objects of two compounds
    Output: mol object of largest common substring
"""
def fmcstimeout(p,q):
    return MCS.FindMCS([p,q], timeout=1).smarts

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

def main():
    # #Test: Time Fragments for 100 strings
    # cpd_smiles = pickle.load(open("../Fragments/100k_RAND_SMILES/100k_smiles_0.p", "rb"))

    ## KEGG smiles (~18k smiles strings, therefore 158 million combinations)
    kegg_smiles = open("Biology/Data/kegg_smiles.txt")
    cpd_smiles = kegg_smiles.readlines()

    cpd_mols = [Chem.MolFromSmiles(smi.strip()) for smi in cpd_smiles]
    cpd_mols = [m for m in cpd_mols if m != None]
    print("Retieved",len(cpd_mols),"random molecules in sample")


    # start = time.time()
    # frags = []
    # for s in (fragments(cpd_mols)): #| loadSmarts(sys.argv[2])):
    #     try:
    #         frags.append((Chem.MolFromSmarts(s),s))
    #     except:
    #         print("AAAGGGHHH",s,"failed")
    # frags=UniqSmarts(frags)
    # print("Found", len(frags), "many fragments")
    # print("Time:", time.time() - start)


    #Define parallel computing
    start = time.time()
    pool = Pool(processes=24)
    cpd_combinations = combinations(cpd_mols, 2)
    frag_smarts = pool.map(findmcs, cpd_combinations)
    frags = []
    for s in (frag_smarts): #| loadSmarts(sys.argv[2])):
        try:
            frags.append((Chem.MolFromSmarts(s),s))
        except:
            pass
    frags=UniqSmarts(frags)
    print("Found", len(frags), "many fragments")
    print("Time:", time.time() - start)

    pickle.dump(frags, open("Biology/Data/KEGG_fragments.p", "wb"))

if __name__ == "__main__":
    main()
