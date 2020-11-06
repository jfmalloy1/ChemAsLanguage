##!/usr/bin/env python
import sys
from rdkit import Chem
from rdkit import RDLogger
from itertools import combinations
from random import random, shuffle, sample
from collections import defaultdict
from rdkit.Chem import MCS
import operator
import datetime
#from timeout import timeout
from math import factorial
from tqdm import tqdm
import json
import pandas as pd
import time
import numpy as np
import matplotlib.pyplot as plt

""" Adopted from Cadaddeu 2014 """

def binomial(n,m):
    return factorial(n) / (factorial(m) * factorial(n-m))

def time_fn(): return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

def findmcs(p,q):
    #@timeout(2)
    def fmcstimeout(p,q):
        return MCS.FindMCS([p,q], timeout=0.1).smarts

    try:
        return fmcstimeout(p,q)
    except:
        print("MCS of", p, "and", q, "timed out.")
        pass

def fragments(mols):
    s=[]
    count = float(binomial(len(mols),2))
    i = 0
    percent = -1
    for p,q in tqdm(combinations(mols,2)):
        try:
            s.append(findmcs(p,q))
        except:
            pass

    return set(s)

def loadSmarts(fn):
    with open(fn,'r') as smartfile:
        return set(smart.strip() for smart in smartfile)


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

def UniqSmarts(frags):
    result = set()
    frags = set((g,t) for (g,t) in frags if g != None) # was successfully converted to smarts
    while frags:
        f,s = frags.pop()
        result.add((f,s))
        frags = set((g,t) for (g,t) in frags
                          if not (sameMolecule((f,s),(g,t)))) # check if they are the same or not. check also on the string to speed up
    return result

"""Creates a dataframe from the Reaxys subset number n between 1-10"""
def read_cpds(n):
    start = time.time()
    df = pd.read_csv("~/Shared/Reaxys/cpd_"+n+".csv", header=None, usecols=[31])

    df.columns=["smiles"]
    end = time.time()
    print("Time for " + n + ":", end-start)
    return df

""" df.apply function for converting smiles to Chem mol object
    Input: smiles string
    Output: Chem mol object, or None if unable
"""
def create_mol(smi):
    try:
        mol = Chem.MolFromSmiles(smi.strip())
    except:
        mol = None
    return mol

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
            if m:
                if m.HasSubstructMatch(f):
                    h[s] += 1
    return h

""" Graphs basic disributions
    Input: h - a dictionary of smarts strings and the number of occurances within KEGG, i: iteration of particular dictionary (e.g., 1-10);
        l - label of graph, rev - True/False distinction for reversability of sorting, fp - filepath for savefig, title - title for plot
    Output: pretty graphs :)
"""
def distribution_graph(h, i, l, rev, fp, title):
    #Calculate AUC to distinguish splits
    yvals = list(sorted(h.values(), reverse=rev))
    xvals = np.linspace(0, 1, num=len(yvals))
    area = np.trapz(yvals, dx=xvals[1])
    plt.plot(xvals, yvals, linewidth=3, label=str(l) + " AUC=" + str(round(area, 2))) #Note: AUC label
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Rank-ordered Compounds")
    plt.ylabel("Occurances")
    plt.title(title)
    plt.legend()


""" RUNTIME NOTES:
    argv[1] = output file
    argv[2] = figure legend label

    argv[1] = random molecues in SMILES format - REMOVED
    argv[2] = more molecules in SMART format (optional?) - REMOVED FOR NOW
"""
if __name__ == "__main__":
    #Remove rdkit warnings
    RDLogger.DisableLog('rdApp.*')

    #Read in full reaction database
    df = pd.DataFrame()
    for i in range(1,11):
        df = df.append(read_cpds(str(i)), ignore_index=True)
        print("Done with subset", i, "...")
        print("Df size", len(df.index))

    #convert all Reaxys to mol object
    start = time.time()
    df = df.sample(frac=0.1) #Sample 1/10th of full database
    df["Mol"] = df["smiles"].apply(create_mol)
    print("Smiles to Mol conversion time:", time.time() - start)

    ## DUPLICATION ##
    for i in range(10): #Note: only doing this once for 1000 cpd test
        print("\n\n------------- ROUND", i, "-------------\n\n")
        #Sample 100 compounds from Reaxys sample
        smiles = df.sample(frac=100/len(df.index))["smiles"].tolist()
        smiles = list(map(str, smiles))

        #Get all of the random molecules as mols
        mols = [Chem.MolFromSmiles(smi.strip()) for smi in smiles]
        mols = [m for m in mols if m != None]
        print("Retieved",len(mols),"random molecules in sample", i)

        #shuffle(mols)
        #Input file 2 - list of molecules in SMART format (still not sure what this adds...)
        frags = []
        #Fragments - determines different fragments within the molecules.
        # NOTE: use a subset of full random mols (mols) or class_mols, as needed
        for s in (fragments(mols)): #| loadSmarts(sys.argv[2])):
            try:
                frags.append((Chem.MolFromSmarts(s),s))
            except:
                print("AAAGGGHHH",s,"failed")
        print("uniquifying")
        #Make sure they are unique
        frags=UniqSmarts(frags)
        print("Found", len(frags), "many fragments")
        h = defaultdict(int)
        count = float(len(mols))

        #Construct histogram over full random molecule set
        print("Constructing histogram of fragments over sample set")
        for m in tqdm(mols):
            for (f,s) in frags:
                if m.HasSubstructMatch(f):
                    h[s] += 1

        print("Writing out fragment occurances V" + str(i))
        #Output file - only smiles fragments
        with open(sys.argv[1] + str(i) + ".txt",'w') as out: #NOTE: CHANGED THIS TO argv[2]
            for k,v in sorted(h.items(), key=operator.itemgetter(1)):
                print(k, file=out)
        #Output file 2 - csv containing smiles fragments & occurances within random molecule set
        with open(sys.argv[1] + str(i) + "_sampleOccurances.csv",'a') as out: #NOTE: CHANGED THIS TO argv[2]
            print("Frags,Occurances", file=out)
            for k,v in sorted(h.items(), key=operator.itemgetter(1)):
                print(str(k) + "," + str(v), file=out)

        ## Find distribution over full Reaxys database ##
        #Find histogram of frag occurances over entire database
        print("Constructing histogram of fragments over full 1/10th of database")
        h = mol_count(df["Mol"].tolist(), frags)

        #Graph things (values, split, label, reverse T/F, plot filepath, plot title)
        distribution_graph(h, i, sys.argv[2], True, sys.argv[2], sys.argv[2].replace("_", " "))
        plt.savefig(sys.argv[1] + "distribution_graph")
        plt.close()
        with open(sys.argv[1] + str(i) + "_fullOccurances.csv",'a') as out: #NOTE: CHANGED THIS TO argv[2]
            print("Frags,Occurances", file=out)
            for k,v in sorted(h.items(), key=operator.itemgetter(1)):
                print(str(k) + "," + str(v), file=out)

        ## Find pre-made distribution over random molecule set ##
        h = pd.read_csv(sys.argv[1] + str(i) + "_sampleOccurances.csv", header=None, skiprows=1, index_col=0, squeeze=True).to_dict()
        distribution_graph(h, 0, sys.argv[2], True, sys.argv[2] + "subset_only", sys.argv[2].replace("_", " ") + " Subset Only")
        plt.savefig(sys.argv[1] + "distribution_graph_subset")
        plt.close()
