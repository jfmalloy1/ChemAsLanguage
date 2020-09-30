from rdkit import Chem
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import numpy as np
from scipy import stats

""" Create a list of mol representations (from rdkit) from a list of smarts strings
    Note: already unique (common_fragments.py takes care of removing duplicates)
    Input: a list of smarts representations
    Output: list of sets: (mol representations, smarts string)
"""
def mol_from_smarts(smarts):
    frags = []
    for s in smarts:
        try:
            frags.append((Chem.MolFromSmarts(s),s))
        except:
            print(s, " failed")

    return frags

""" Creates a histogram of occurances of fragments within the overall cpd distribution
    Input:  mols - mol representations of overall distribution
            frags - fragments found within KEGG
    Output: dictionary containing number of time a fragment appears in the overall distribution
"""
def mol_count(mols, frags):
    h = {}
    for (f,s) in tqdm(frags):
        h[s] = 0
        for m in mols:
            if m.HasSubstructMatch(f):
                h[s] += 1

    return h

""" Find the base fragments - those that appear in all splits
    Input: list of lists of fragments. Overarching list - all splits, each sublist - fragments within that split
    Output: fragments which are common to all fragment splits
"""
def base_frags(frags):
    # #Initial test: how many fragments do split 0 & 1 have?
    # set0 = set(frags[0])
    # print(frags[0])
    # set1 = set(frags[1])
    # print("\n\n\n")
    # print(frags[1])
    # base_frags = set0.intersection(set1)
    # print(base_frags)
    # print(len(base_frags))

    #Initial step - base fragments are the initial split
    base_frags = frags[0]
    #Total number of fragments - starting with the initial fragment set
    total_frags = frags[0]

    #Statistics over number of frags
    frag_stats = [len(frags[0])]

    #Find the intersection of all further splits
    for i in range(1, len(frags)):
        print("Length of split", i, "is:", len(frags[i]))
        frag_stats.append(len(frags[i]))
        base_frags = set(frags[i]).intersection(set(base_frags))
        total_frags = list(set(total_frags + frags[i]))

    print("Number of base fragments:", len(base_frags))
    # print(base_frags)

    #Goal - find total number of different fragments
    print("Number of total fragments:", len(total_frags))

    print("Fragment mean:", np.mean(frag_stats))
    print("Fragment std:", np.std(frag_stats))

""" Graphs basic disributions
    Input: h - a dictionary of smarts strings and the number of occurances within KEGG
    Output: pretty graphs :)
"""
def distribution_graph(h, i):
    #Calculate AUC to distinguish splits
    yvals = list(sorted(h.values(), reverse=True))
    xvals = np.linspace(0, 1, num=len(yvals))
    area = np.trapz(yvals, dx=xvals[1])
    plt.plot(xvals, yvals, label = "Split " + str(i) + " AUC=" + str(round(area, 2)))
    plt.legend()

def main():
    ## Read in mol files of KEGG/domain/LUCA ##
    with open("../kegg_smiles.txt",'r') as smiles:
        mols = [Chem.MolFromSmiles(smi.strip()) for smi in smiles]
        mols = [m for m in mols if m != None]

    #Read in pre-generated fragments from a file
    frags = []
    dirpath = "Tests/Twenty_cpds_nomix/"
    for file in os.listdir(dirpath):
        fp = dirpath + file
        with open(fp) as f:
            frags.append([line.rstrip('\n') for line in f])

    # ## Find repeatability ##
    # base_frags(frags)

    ## Find distribution ##
    for i in range(len(frags)):
        #Turn smarts list into list of sets (mol, smarts)
        f = mol_from_smarts(frags[i])
        #Find histogram of frag occurances
        h = mol_count(mols, f)
        #Graph things
        distribution_graph(h, i)

    plt.show()


if __name__ == "__main__":
    main()
