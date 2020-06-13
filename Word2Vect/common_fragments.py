##!/usr/bin/env python
import sys
from rdkit import Chem
from itertools import combinations
from random import random, shuffle
from collections import defaultdict
from rdkit.Chem import MCS
import operator
import datetime
#from timeout import timeout
from math import factorial

""" Adopted from Cadaddeu 2014 """

def binomial(n,m):
    return factorial(n) / (factorial(m) * factorial(n-m))

def time(): return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

def findmcs(p,q):
    #@timeout(2)
    def fmcstimeout(p,q):
        return MCS.FindMCS([p,q]).smarts

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
    for p,q in combinations(mols,2):
        try:
            s.append(findmcs(p,q))
        except:
            pass

        p = int(i/count * 100)
        if p > percent:
            percent = p
            print(percent, time())
        i += 1
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

if __name__ == "__main__":
    #Input file - list of random molecules in SMILES format
    print("Getting random molecules from", sys.argv[1])
    with open(sys.argv[1],'r') as smiles:
        mols = [Chem.MolFromSmiles(smi.strip()) for smi in smiles]
        mols = [m for m in mols if m != None]
        shuffle(mols)
    print("Retieved",len(mols),"random molecules")
    print("Aquiring random molecule fragments (and combining with molecules from", sys.argv[2])

    #Input file 2 - list of molecules in SMART format (still not sure what this adds...)
    frags = []
    #Fragments - determines different fragments within the molecules. Why are only the first 20 random molecules used?
    for s in (fragments(mols[:100]) | loadSmarts(sys.argv[2])):
        try:
            frags.append((Chem.MolFromSmarts(s),s))
        except:
            print("AAAGGGHHH",s,"failed")
    print("uniquifying")
    #Make sure they are unique
    frags=UniqSmarts(frags)
    print("Found", len(frags), "many fragments")
    h = defaultdict(int)
    percent = -1
    i = 0
    count = float(len(mols))

    print("Constructing histogram of fragments")
    for m in mols:
        for (f,s) in frags:
            if m.HasSubstructMatch(f):
                h[s] += 1
        p = int(i / count * 100)
        if p > percent:
            percent = p
            print(percent, time())
        i += 1

    print("Writing out histogram")
    #Output file
    with open(sys.argv[3],'w') as out:
        for k,v in sorted(h.items(), key=operator.itemgetter(1)):
            print(k, file=out)
