import pandas as pd

def main():
    #Read in smiles strings
    fp = "../../ChiralityXP/kegg_chirality.csv"
    df = pd.read_csv(fp)

    #Write smiles to a text file
    for item, row in df.iterrows():
        print(row["Original SMILES"], file=open("kegg_smiles.txt", "a"))

if __name__ == "__main__":
    main()
