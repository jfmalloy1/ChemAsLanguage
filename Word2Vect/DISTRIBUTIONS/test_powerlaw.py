import powerlaw
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from rdkit import Chem
import pandas as pd
from tqdm import tqdm
import os
import sys
import math

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

""" Graphs basic disributions
    Input: h - a dictionary of smarts strings and the number of occurances within KEGG
    Output: pretty graphs :)
"""
def distribution_graph(h, i):
    #Calculate AUC to distinguish splits
    yvals = list(sorted(h.values(), reverse=True))
    #xvals = np.linspace(0, 1, num=len(yvals)) #NOTE: for AUC normalization calculations
    #area = np.trapz(yvals, dx=xvals[1]) #NOTE: for AUC calculations
    xvals = np.arange(len(yvals))
    plt.plot(xvals, yvals, label = "Split " + str(i)) #NOTE:for AUC + " AUC=" + str(round(area, 2)))
    # plt.xscale("log")
    # plt.yscale("log")
    plt.legend()

#NOTE: sys.argv[0] - name of directory within Tests to analyze
def main():
    # ## Read in mol files of KEGG/domain/LUCA ##
    # with open("../kegg_smiles.txt",'r') as smiles:
    #     mols = [Chem.MolFromSmiles(smi.strip()) for smi in smiles]
    #     mols = [m for m in mols if m != None]
    #
    # #Read in pre-generated fragments from a file
    # label = sys.argv[1]
    # frags = []
    # fp = "Tests/" + label + "/dup_100cpds_0.txt"
    # with open(fp) as f:
    #     frags.append([line.rstrip('\n') for line in f])
    #
    # ## Find distribution ##
    # for i in range(len(frags)):
    #     #Turn smarts list into list of sets (mol, smarts)
    #     f = mol_from_smarts(frags[i])
    #     #Find histogram of frag occurances
    #     h = mol_count(mols, f)
    #     # #Graph things
    #     # distribution_graph(h, i)

    ## Find pre-made distribution over random molecule set ##
    h = pd.read_csv("Tests/Hundred_cpds_random_subsampleOccurances/dup_0_occurances.csv", header=None, skiprows=1, index_col=0, squeeze=True).to_dict()
    try:
        del(h["Frags"])
    except:
        pass

    #powerlaw fit
    values = list(map(float, list(h.values())))
    #print(values)
    fit = powerlaw.Fit(values, xmin=1)
    #print(label)
    print("alpha:", fit.alpha)
    print("Standard Error:", fit.sigma)
    a = float(fit.alpha)
    print("KS Distance", fit.D)

    print(fit.sigma/fit.alpha)
    bin_edges, probability = fit.pdf()
    print(probability)
    print(len(probability))
    print(len(values))

    #Convert D to alpha (??)
    p = 2*math.exp((fit.D**2 * len(values) * 2 * len(probability)) / (1 + (len(probability) / len(values))))
    print(p)

    # #Testing - generate rvs from a powerlaw distribution with the fitted alpha
    # pl_rvs = stats.powerlaw.rvs(fit.alpha, size=len(values))
    #
    # print(stats.kstest(probability, 'powerlaw', args=(fit.alpha, 1)))
    # # print(pl_rvs)
    # # plt.plot(sorted(probability), label="modeled probability")
    # # plt.plot(sorted(pl_rvs), label="Fitted rvs")
    # # plt.legend()
    # # plt.show()
    # print("Stats KS:", stats.ks_2samp(values, pl_rvs))
    # # plt.plot(bin_edges[1:], probability)
    # # plt.plot(sorted(h.values(), reverse=True))
    # # plt.show()


    # #Plotting: taken from https://stackoverflow.com/questions/16640496/python-plot-and-powerlaw-fit
    # figPDF = fit.power_law.plot_pdf(color="darkgreen", linestyle="--", label="Power Law PDF\nKS distance= " + str(round(fit.D,2)))
    # fit.plot_pdf(color="darkgreen", linewidth=2, label="Fragment Distribution pdf", ax=figPDF)
    # fit.power_law.plot_ccdf(color="slateblue", linestyle="--", label="Power Law CCDF", ax=figPDF)
    # fit.plot_ccdf(color="slateblue", linewidth=2, label="Fragment Distribution CCDF", ax=figPDF)
    #
    # plt.legend()
    # plt.xlabel("Fragment Occurance")
    # plt.ylabel("Probability (occurance)")
    # #plt.savefig("Annual_Review_PLTest")
    # plt.show()


    # data = stats.powerlaw.rvs(3, size=1000)
    #
    # fit = powerlaw.Fit(data)
    # print(fit.alpha)
    # print(fit.D)
    #
    # powerlaw.plot_pdf(data, color="b")
    # plt.plot(data)
    # plt.show()
    # #plt.savefig("testPwrlaw")
    #
    # print(fit.distribution_compare("power_law", "exponential", normalized_ratio=True))
    #
    # # plt.hist(data, bins=20)
    # # pwrlaw = stats.powerlaw(fit.power_law.alpha)
    # # x = np.linspace(0,1,100)
    # # plt.plot(x, pwrlaw.pdf(x))
    # # plt.show()

if __name__ == "__main__":
    main()
