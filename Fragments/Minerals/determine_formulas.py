import pandas as pd
import chemparse
import re

def determine_formulas():
    pass

def main():
    df = pd.read_csv("Data/IMA_abiotic_labels.csv")
    formulas = df["IMA Chemistry (plain)"].tolist()

    dash_count = 0
    extra_count = 0
    for f in formulas:
        if '+' in f:
            f = re.sub('\d\+','',f)
            #print(chemparse.parse_formula(f))
        elif '-' in f:
            print(f)
            dash_count += 1
        elif re.findall('[^\da-zA-Z()\\box]+', f):
            print(f)
            print(chemparse.parse_formula(f))
            extra_count += 1
    print("Dash count:", dash_count)
    print("Extra count:", extra_count)

if __name__ == "__main__":
    main()
