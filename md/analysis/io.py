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
