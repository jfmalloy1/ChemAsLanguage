import pandas as pd
import numpy as np

def main():
    #read in V1 coocurrance
    df = pd.read_csv("frags_coocurrance_V1.csv")

    #convert df to numpy array (drop first column, convert)
    sub_df = df.drop(df.columns[0], axis=1)
    X = sub_df.to_numpy()
    print(X)

    #transpose matrix
    X = X + X.T - np.diag(np.diag(X))
    #print(X)

    #Write back out to csv file
    df_transposed = pd.DataFrame(data = X, index = df.columns[1:], columns = df.columns[1:])
    #print(df_transposed)
    df_transposed.to_csv("frags_coocurrance_full.csv")


if __name__ == "__main__":
    main()
