import pandas as pd
from tqdm import tqdm

""" Filters the IMA database by age
    Input: IMA dataframe, age to filter by (in MYA ago), filename of new file
    Output: file (written to the Data/ directory) containing the filtered dataframe
"""
def filter_by_age(df, age, filename):
    sub_df = df[df["Oldest Known Age (Ma)"] >= age]
    print("Number of minerals older than", age, ":", len(sub_df))
    sub_df.to_csv("Data/"+filename+".csv", index=False)

def main():
    #Read in all minerals from IMA
    df = pd.read_csv("Data/IMA_ALL.csv")
    print("Length of full dataframe: ", len(df))

    filter_by_age(df, 3500, "IMA_3500mya")

if __name__ == "__main__":
    main()
