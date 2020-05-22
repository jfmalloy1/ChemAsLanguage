import pandas as pd
import pubchempy as pcp
import operator

#Return the cas number and the name of opioids
def get_info(path):
    df = pd.read_csv(path)
    return list(df["cas_number"]), list(df["name"])

#Returns the pubchem ID of a compound (based on name only)
#Only returns the first cpd_id taken from Pubchem search
def get_pubchem_id(name):
    try:
        cpd_id = pcp.get_cids(name, "name")
        return cpd_id[0]
    except:
        return ""

def get_molecular_weight(name):
    try:
        cpd = pcp.get_compounds(name, "name")
        return cpd[0].to_dict(properties=["molecular_weight"])["molecular_weight"]
    except:
        return 300 #not real mw, but it's in the middle - therefore won't be misued in min/max calculations

def main():
    #Read in drugbank csv, return cas_ids
    cas_ids, names = get_info("../CitationNetworks/drugs/opioids_filtered.csv")
    #Find all pubchem ids of opioids
    opioids_pubchem_ids = []
    opioids_mw = {}
    for name in names:
        print(name)
        #Find all pubchem ids - needed for citation networks
        opioids_pubchem_ids.append(get_pubchem_id(name))
        #Find all molecular weights
        opioids_mw[name] = get_molecular_weight(name)

    #Find all molecular weights of opioids
    #Goal = find smallest and largest
    print(opioids_mw)
    print(max(opioids_mw.keys(), key=lambda k:opioids_mw[k]))
    print(min(opioids_mw.keys(), key=lambda k:opioids_mw[k]))

if __name__ == "__main__":
    main()
