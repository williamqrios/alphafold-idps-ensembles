from pathlib import Path

import numpy as np
import pandas as pd


def load_rg(path: Path) -> pd.DataFrame:
    """
    Loads the radius of gyration data from a CSV file and returns it as a DataFrame.
    """
    if not path.exists():
        raise FileNotFoundError(f"Radius of gyration file not found: {path}")

    df = pd.read_csv(path)
    return df


def load_dssp(
    path: Path,
    time_array: np.ndarray | None = None,
) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"DSSP file not found: {path}")

    dssp_data = np.load(path, allow_pickle=True)
    df = pd.DataFrame(
        dssp_data, columns=[f"residue_{i + 1}" for i in range(dssp_data.shape[1])]
    )
    if time_array is not None:
        df.insert(0, "time", time_array)

    return df


def process_dssp(
    dssp_df: pd.DataFrame, code: str = "H", timeseries: bool = True
) -> pd.DataFrame:
    """
    code: DSSP code, 'H' for helix, 'C' for coil, 'B' for beta sheets.
    Option "timeseries" calculates the fraction of the structure per timeframe to plot it as a time series.
    Setting "timeseries" to False returns global averages per residue, per conformation.
    """
    if code not in ["H", "C", "E"]:
        raise ValueError(
            "Invalid DSSP code. Use 'H' for helix, 'C' for coil, or 'E' for beta sheets."
        )
    axis = 1 if timeseries else 0

    structure_series: pd.Series = dssp_df.iloc[:, 1:].apply(
        lambda row: (row == code).astype(int).mean(), axis=axis
    )  # type: ignore
    # Rebuild the df
    structure_df: pd.DataFrame = pd.DataFrame(
        {
            "time" if timeseries else "residue": dssp_df["time"]
            if timeseries
            else structure_series.index + 1,
            "fraction" if timeseries else "average": structure_series.values,
        }
    )
    return structure_df


def load_exp_chemical_shifts(path: Path) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Assumes CSV file with columns: res, atom_name, value
    """
    df = pd.read_csv(path)
    available_shifts = df.atom_name.unique()

    return df, available_shifts


def load_chemical_shifts(
    path: Path, available_shifts: np.ndarray, fraction_discarded: float = 0.0
) -> pd.DataFrame:
    """Loads estimated chemical shifts, filtering out atoms with no reported experimental data."""

    # Use only chemical shifts that are available in experimental data
    # available_shifts = self.cs_exp.atom_name.unique()
    if "HA" in available_shifts:
        # HA2 corresponds to glycine's other hydrogen alpha
        available_shifts = np.append(available_shifts, "HA2")

    cs_data = pd.read_pickle(path)
    # Discard partion of the trajectory
    mask = (cs_data.atom_name.isin(available_shifts)) & (
        cs_data.time > cs_data.time.max() * fraction_discarded
    )
    cs_data = cs_data[mask].reset_index(drop=True).replace({"atom_name": {"HA2": "HA"}})
    return cs_data



def get_cs_per_frame(df_sim: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """
    Formats dataframe of chemical shifts into a dictionary of dataframes, one for each atom type, with time as index and residues as columns.
    """
    output = {}
    shift_types = df_sim.atom_name.unique()
    for shift_type in shift_types:
        df_sim_shift = (
            df_sim[df_sim.atom_name == shift_type]
            .drop(columns=["atom_name"])
            .pivot_table(index="time", columns="res", values="value")
            .reset_index()
        )
        output[shift_type] = df_sim_shift
    return output
