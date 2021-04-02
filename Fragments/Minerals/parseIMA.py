import pandas as pd
from tqdm import tqdm
import numpy as np
import re

""" Parse the .js file (saved as .txt) of IMA minerals
    Input: filepath to IMA text file
    Output: csv file (saved in Data folder) of relevant data
"""
def parse_data(fp):
    #Initalize pandas dataframe
    df = pd.DataFrame(columns=["mineral name", "official formula", "chemistry", "a", "b", "c", "alpha", "beta",
        "gamma"])
    #Read the 1st 5 lines of IMA.txt
    with open(fp) as f:
        next(f)
        for i in tqdm(range(22384)): #There are 22385 lines in IMA.txt, so skipping the first one
            line = f.readline().split('|')
            l = len(df)
            #df = df.append(pd.Series, index = )
            df.loc[l] = line[2:11]
            #print(len(line.split('|')))

    print(df)
    df.to_csv("Data/IMA.csv")

""" Removes the parentheses, as well as all characters within those parentheses, of a string
    Necessary to parse a/b/c columns in IMA dataframe
    Input: string, regex
    Output: string with no parentheses
"""
def remove_parentheses(s, regex):
    if isinstance(s, str):
        if '(' in s:
            s = re.sub(regex, "", s)
            return s
        else:
            return s
    else:
        return s

""" Remove error tolerances (represented by (error)) in a/b/c coordinates
    Input: file path of IMA csv
    Output: IMA csv with error tolerances removed
"""
def remove_error_tol(fp):
    df = pd.read_csv(fp)
    regex = re.compile("[\(\[].*?[\)\]]") #from https://stackoverflow.com/questions/14596884/remove-text-between-and
    for index, row in tqdm(df.iterrows()):
        for col in ["a", "b", "c"]:
            df.loc[index, col] = remove_parentheses(row[col], regex)

    return df

""" Find minerals (within a given tolerance) given a/b/c coordinates
    Input: filepath to IMA csv, tolerance percentage (in 0-1 limits), a, b, c coordinates
    Output: list of minerals which fit the tolerance parameters
"""
def find_abc(fp, tol, a_val, b_val, c_val):
    df = pd.read_csv(fp)
    print(len(df))

    #Change a/b/c to numerics
    df["a"] = pd.to_numeric(df["a"], errors="coerce")
    df["b"] = pd.to_numeric(df["b"], errors="coerce")
    df["c"] = pd.to_numeric(df["c"], errors="coerce")

    #Error tolerance from https://stackoverflow.com/questions/41686930/search-for-value-within-a-range-in-a-pandas-dataframe
    sub_df = df[df[["a"]].apply(np.isclose, b=a_val, atol=0.1).any(1)]
    sub_df = sub_df[sub_df[["b"]].apply(np.isclose, b=b_val, atol=0.1).any(1)]
    sub_df = sub_df[sub_df[["c"]].apply(np.isclose, b=c_val, atol=0.1).any(1)]

    print(len(sub_df))
    print(sub_df)


def main():
    # #Parse the IMA database
    # parse_data("Data/IMA.txt")
    # #Remove error tolerances within a/b/c data
    # df = remove_error_tol("Data/IMA.csv")
    # df.to_csv("Data/IMA_cleaned.csv")

    #Find a mineral given a/b/c - using Abellaite as a test
    find_abc("Data/IMA_cleaned.csv", 0.1, 5.26, 5.26, 13.463)

if __name__ == "__main__":
    main()
