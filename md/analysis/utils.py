import os

import pandas as pd


def unpivot_df(df):
    """
    Unpivot a DataFrame from wide to long format.
    """
    df_unpivoted = df.reset_index().melt(
        id_vars=["resSeq", "name"], var_name="time", value_name="value"
    )
    df_unpivoted.rename({"resSeq": "res", "name": "atom_name"}, axis=1, inplace=True)
    return df_unpivoted


def process_shift_chunks(name: str, conformation: str, forcefield: str):
    """
    Merges multiple pickle files containing chemical shift data into a single DataFrame and saves it as a new pickle file.
    """
    path = f"{name}/{conformation}/{forcefield}"
    pkl_files = sorted([f for f in os.listdir(path) if f.endswith(".pkl")])
    merged = pd.DataFrame()
    for pkl in pkl_files:
        data = unpivot_df(pd.read_pickle(f"{path}/{pkl}"))
        if merged.empty:
            merged = data
        else:
            merged = pd.concat([merged, data], ignore_index=True)
    merged.to_pickle(f"{path}/shifts_merged.pkl")


# Example usage:
def main():
    protein_name = "p53n"
    conformation = "extended"
    forcefield = "a99disp"
    process_shift_chunks(
        name=protein_name, conformation=conformation, forcefield=forcefield
    )


if __name__ == "__main__":
    main()
