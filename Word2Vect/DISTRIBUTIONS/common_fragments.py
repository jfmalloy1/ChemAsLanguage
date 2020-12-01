##!/usr/bin/env python
import sys
from rdkit import Chem
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
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import numpy as np
from scipy import stats

""" Adopted from Cadaddeu 2014 """

def binomial(n,m):
    return factorial(n) / (factorial(m) * factorial(n-m))

def time(): return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

def findmcs(p,q, t):
    #@timeout(2)
    def fmcstimeout(p,q):
        return MCS.FindMCS([p,q], timeout=1).smarts

    try:
        return fmcstimeout(p,q)
    except:
        print("MCS of", p, "and", q, "timed out.")
        pass

def fragments(mols, t):
    s=[]
    count = float(binomial(len(mols),2))
    i = 0
    percent = -1
    for p,q in tqdm(combinations(mols,2)):
        try:
            s.append(findmcs(p,q, t))
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

""" Find specific compound classes given the KEGG Brite json file, given a specific labl
    Input: filepath to the KEGG BRITE classification, KEGG compound list (including smiles), classification label to find
    Output: SMILES strings associated with each compound class
"""
def cpd_classes(brite_fp, cpd_fp, label):
    with open(brite_fp) as json_file:
        br_json = json.load(json_file)

    #Create a dictionary with each compound having the appropriate label
    classes = []
    cpd_class = {}
    for key1 in br_json["children"]:
      classes.append(key1["name"])
      for key2 in key1["children"]:
        for key3 in key2["children"]:
          for cpd in key3["children"]:
            cpd_class[cpd["name"][:6]] = key1["name"]
            #print(cpd["name"][:6] + " " + key1["name"])

    #Read in all KEGG compound labels
    df = pd.read_csv(cpd_fp)

    #Find & return all smiles strings associated with a specific label
    cpds = [k for k, v in cpd_class.items() if v == label]
    return df[df["C"].isin(cpds)]["S"].tolist()

""" Read in a specific percentile of universal compounds (5 - top 5%, 10 - top 10%, etc...)
    Input: percentile of compound (must be 5, 10, or 15)
    Output: list of smiles corresponding to unviversal compounds
"""
def percentiles(p):
    #Read in compound labels
    cpds = []
    if p == 5:
        with open("95percentile.txt") as f:
            cpds = [line.rstrip("\n") for line in f]
    elif p == 10:
        with open("90percentile.txt") as f:
            cpds = [line.rstrip("\n") for line in f]
    else:
        with open("85percentile.txt") as f:
            cpds = [line.rstrip("\n") for line in f]

    #Read in all KEGG compound smiles data
    df = pd.read_csv("../chiral_molweight_formula_labels.csv")
    #Return list of smiles
    smiles = list(map(str, df[df["C"].isin(cpds)]["S"].tolist()))
    return smiles

""" Create a list of mol representations (from rdkit) from a list of smarts strings
    Note: already unique (common_fragments.py takes care of removing duplicates)
    Input: a list of (mol, smarts) representations
    Output: list of sets: (mol representations, smarts string)
"""
def mol_from_smarts(smarts):
    mols = []
    for m, s in smarts:
        mols.append(m)

    return mols

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

""" Find the base fragments - those that appear in all splits
    Input: list of lists of fragments. Overarching list - all splits, each sublist - fragments within that split
    Output: fragments which are common to all fragment splits
"""
def base_frags(frags):
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

""" Find all smiles strings associated with a particular domain
    Input: Domain from ["Eukarya", "Bacteria", "Archaea"]
    Output: list of all smiles strings associated with that particular domain
"""
def domain_smiles(domain):
    fp = ""
    if domain == "Eukarya":
        fp = "eukarya_cpds.csv"
    elif domain == "Bacteria":
        fp = "bacteria_cpds.csv"
    elif domain == "Archaea":
        fp = "archaea_cpds.csv"

    domain_df = pd.read_csv(fp)
    #kegg_data_curated is the same as chiral_molweight_formula_labels, containing smiles in the "S" column
    kegg_df = pd.read_csv("kegg_data_curated.csv") #Assumes kegg_data_curated.csv is in above directory

    #Return smiles strings (in list object) of that particular domain
    return kegg_df[kegg_df["C"].isin(domain_df["compounds"])].dropna(subset=["S"])["S"].tolist()

""" Find all smiles compounds associated with minerals (COD database subset)
    Input: filepath to csv file containing a subset of the COD database
    Output: List of smiles
"""
def mineral_smiles(fp):
    df = pd.read_csv(fp)
    return df["SMI"].tolist()

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
    argv[1] = random molecues in SMILES format
    argv[2] = output filepath (ex: "Tests/Hundred_cpds/")
    argv[3] = figure label

    argv[2] = more molecules in SMART format (optional?) - REMOVED FOR NOW
"""
if __name__ == "__main__":
    #Input file - list of random molecules in SMILES format
    print("Getting (database) compounds from", sys.argv[1])

    #Make directory (if it does not already exist)
    if not os.path.isdir(sys.argv[2]):
        os.mkdir(sys.argv[2])

    # ## Compound Classes ##
    # smiles = cpd_classes("../br08001.json", "../chiral_molweight_formula_labels.csv", "Hormones and transmitters")
    # mols = [Chem.MolFromSmiles(smi.strip()) for smi in smiles]
    # mols = [m for m in mols if m != None]
    # print("Retieved",len(mols),"random molecules")

    # ## Percentiles ##
    # #Input should be 5, 10, or 15 to refer to top 5%, 10%, and 15% of universal compounds
    # smiles = percentiles(5)
    # mols = [Chem.MolFromSmiles(smi.strip()) for smi in smiles]
    # mols = [m for m in mols if m != None]
    # print("Retieved",len(mols),"classified compounds")

    # ## Domains ##
    # #Goal - have either Eukarya, Bacteria, or Archaea compounds (general compounds) be used for fragments
    # smiles = domain_smiles("Eukarya")
    # mols = [Chem.MolFromSmiles(smi.strip()) for smi in smiles]
    # mols = [m for m in mols if m != None]
    # print("Retieved",len(mols),"random molecules")

    ## Minerals
    #Goal - see mineral occurance, both within minerals themselves and across KEGG (for now)
    smiles = mineral_smiles("COD_SMI_IMA_subset218.csv")
    mols = [Chem.MolFromSmiles(smi.strip()) for smi in smiles]
    mols = [m for m in mols if m != None]
    print("Retieved",len(mols),"random molecules")

    #Get all of KEGG
    with open(sys.argv[1],'r') as db_smiles:
        db_mols = [Chem.MolFromSmiles(smi.strip()) for smi in db_smiles]
        db_mols = [m for m in db_mols if m != None]
    print("Retieved",len(db_mols),"base molecules")
    #print("Aquiring random molecule fragments") # (and combining with molecules from", sys.argv[2])

    ## DUPLICATION ##
    for i in range(1):
        #shuffle(mols)
        #Input file 2 - list of molecules in SMART format (still not sure what this adds...)
        frags = []
        #Fragments - determines different fragments within the molecules.
        # # NOTE: use a subset of full random mols (mols) or class_mols, as needed
        #mols = sample(db, 100)
        #mols = db_mols #Note: for full database test
        # for t in [0.01, 0.1, 1, 10, 100]: #Note: testing timeout for MCS algorithm
        t = 0.1 #timeout time for MCS - 0.1
        for s in (fragments(mols, t)): # either mols or a sample, depending on if a subsample is taken or not
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
        print("Constructing histogram of fragments")
        for m in tqdm(mols): #mol_sample or mols, depending on if a subsample is taken or not
            for (f,s) in frags:
                if m.HasSubstructMatch(f):
                    h[s] += 1

        print("Writing out histogram V" + str(i))
        #Output file - only smiles fragments
        with open(sys.argv[2] + str(i) + ".txt",'w') as out: #NOTE: CHANGED THIS TO argv[2]
            for k,v in sorted(h.items(), key=operator.itemgetter(1)):
                print(k, file=out)
        #Output file 2 - csv containing smiles fragments & occurances within random molecule set
        with open(sys.argv[2] + str(i) + "_sampleOccurances.csv",'a') as out: #NOTE: CHANGED THIS TO argv[2]
            print("Frags,Occurances", file=out)
            for k,v in sorted(h.items(), key=operator.itemgetter(1)):
                print(str(k) + "," + str(v), file=out)

        ## DISTRIBUTION ANALYSIS ##

        # ## Find repeatability ## #Note: if more than one loop is done, this (and everything below) should be outside
        # base_frags(frags)

        ## Find distribution over full database ##
        #Find histogram of frag occurances over entire database
        h = mol_count(db_mols, frags)

        #Graph things (values, split, label, reverse T/F, plot filepath, plot title)
        distribution_graph(h, i, sys.argv[3], True, sys.argv[2], sys.argv[3].replace("_", " "))
        with open(sys.argv[2] + str(i) + "_fullOccurances.csv",'a') as out: #NOTE: CHANGED THIS TO argv[2]
            print("Frags,Occurances", file=out)
            for k,v in sorted(h.items(), key=operator.itemgetter(1)):
                print(str(k) + "," + str(v), file=out)

    ## Find pre-made distribution over random molecule set ##
    h = pd.read_csv(sys.argv[2] + str(i) + "_fullOccurances.csv", header=None, skiprows=1, index_col=0, squeeze=True).to_dict()
    distribution_graph(h, 0, sys.argv[3], True, sys.argv[2] + "subset_only", sys.argv[3].replace("_", " ") + " Subset Only")

    plt.savefig(sys.argv[2] + "distribution_graph")
    plt.close()
