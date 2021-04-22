from rdkit import Chem
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import numpy as np
from scipy import stats
import pickle
import random
from multiprocessing import Pool
import time

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
            pass

    return frags

""" Creates a histogram of occurances of fragments within the overall cpd distribution
    Input:  mols - mol representations of overall distribution
            frags - fragments found within KEGG
    Output: dictionary containing number of time a fragment appears in the overall distribution {smarts:occurances}
"""
def mol_count(mols, frags):
    h = {}
    for (f,s) in tqdm(frags):
        h[s] = 0
        for m in mols:
            if m.HasSubstructMatch(f):
                h[s] += 1
    return h

""" Parallel version of mol_count
    Inputs: list of mol representations of overall distribution, individual fragment to be counted
    Output: list of smarts string & count (will be converted to dictionary later)
"""
def mol_count_parallel(mols, frag):
    count = 0
    f,s = frag
    for m in mols:
        if m.HasSubstructMatch(f):
            count += 1

    return [s, count]

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
    Input: h - a dictionary of smarts strings and the number of occurances within KEGG, i: iteration of particular dictionary (e.g., 1-10)
    Output: pretty graphs :)
"""
def distribution_graph(h, i):
    #Calculate AUC to distinguish splits
    yvals = list(sorted(h.values(), reverse=True))
    xvals = np.linspace(0, 1, num=len(yvals))
    #area = np.trapz(yvals, dx=xvals[1])
    plt.plot(xvals, yvals, label = "100 Random KEGG Compounds", color="darkgreen", linewidth=3)#"Split " + str(i) + " AUC=" + str(round(area, 2))) #Note: AUC label
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Rank-ordered Compounds")
    plt.ylabel("Occurances in KEGG")
    plt.legend
    plt.show()

""" Convert a pickled distribution to csv format
    Input: File path to a pickled dictionary
    Output: csv file (same location) of data
"""
def convert_dist_toCSV(fp):
    # for label in ["5000cpds"]:#["1000cpds", "2000cpds", "3000cpds", "4000cpds"]:#, "5000cpds"]:
    #     print(label)
    sample_frags = pickle.load(open(fp, "rb"))
    df = pd.DataFrame.from_dict(sample_frags, orient="index", columns=["Occurances"])
    df["Fragments"] = df.index
    df = df.reset_index(drop=True)
    df.to_csv(fp[:-2] + ".csv")

def main():
    # # ## Read in mol objects of KEGG ##
    # with open("Biology/Data/kegg_smiles.txt",'r') as smiles:
    #     mols = [Chem.MolFromSmiles(smi.strip()) for smi in smiles]
    #     mols = [m for m in mols if m != None]
    #
    # # Parallel Occurances calculations
    # pool = Pool(processes=10)
    # dirpath = "Biology/Data/Tests/Timeout/"
    # for file in os.listdir(dirpath): #For reading in all fragments
    #     if file.endswith("_unique.p"): #ensure only unique fragment sets are counted
    #         start = time.time()
    #         fp = dirpath + file
    #         print("Analyzing:", fp)
    #         frags = pickle.load(open(fp, "rb")) #Load in fragments
    #         print("Analyzing", len(frags), "fragments")
    #         frag_occurances = []
    #         frag_occurances = pool.starmap(mol_count_parallel, tqdm(zip([mols]*len(frags), frags), total=len(frags))) #Parallel occurances calculations
    #         h = {f[0]: f[1] for f in frag_occurances} #Convert list of lists into dictionary
    #         pickle.dump(h, open(dirpath + file[:-2] + "_occurances.p", "wb")) #Save dictionary to new pickle file
    #         print("Time:", time.time() - start)
    #         print()

    # with open("Tests/Hundred_cpds/dup_100cpds_0.txt") as f:
    #     frags.append([line.rstrip('\n') for line in f])
    #
    # # ## Find repeatability ##
    # # base_frags(frags)

    # ## Find distribution of a fragment sample over full database ##
    # for label in ["1000cpds", "2000cpds", "3000cpds", "4000cpds", "5000cpds"]:
    #     print(label)
    # frags = pickle.load(open("Biology/Data/KEGG_fragments_full.p", "rb"))
    # h = mol_count(mols, frags)
    #
    # pickle.dump(h, open("Biology/Data/KEGG_fragments_full_occurances.p", "wb"))

    for f in os.listdir("Biology/Data/Tests/Timeout/"):
        if f.endswith("_occurances.p"):
            convert_dist_toCSV("Biology/Data/Tests/Timeout/" + f)

    # ## Find pre-made distribution over random molecule set ##
    # h = pd.read_csv("Tests/Hundred_cpds_random_subsampleOccurances/dup_0_occurances.csv", header=None, skiprows=1, index_col=0, squeeze=True).to_dict()
    # print(h)
    # distribution_graph(h, 0)
    #
    # plt.show()


if __name__ == "__main__":
    main()
