from tqdm import tqdm

def main():
    #Open mineral database
    for line in tqdm(open("all_cod_smiles.txt")):
        print(line.split("\t")[0], file=open("mineral_smiles.txt", "a"))

if __name__ == "__main__":
    main()
