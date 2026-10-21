from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd
import yaml
from scipy import optimize

MU0 = 4.0e-7 * np.pi
HBAR = 1.055e-34
GAMMA_H = 267.513e6  # rad s^-1 T^-1
GAMMA_N = -27.166e6  # rad s^-1 T^-1
NH_DISTANCE = 0.101e-9  # m
CSA_N = -160.0e-6


def load_t0() -> dict[str, dict[str, dict[str, float]]]:
    """Load equilibration time data from a YAML file."""
    config_path = Path(__file__).parent.joinpath("equilibration_time.yaml")
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def load_traj(
    protein_name: str,
    conformation: str,
    forcefield: str,
    traj_name: str = "md_protein.xtc",
    top_name: str = "md_protein.pdb",
    t0: float | None = None,
) -> md.Trajectory:
    """
    Load a trajectory using MDTraj.
    This includes only the protein, already centered in the box and free of PBC effects.
    Each frame is separated from the next frame by time step of 200 ps (0.2 ns).
    """
    traj_path = Path(__file__).parent.joinpath(
        protein_name, conformation, forcefield, traj_name
    )
    top_path = Path(__file__).parent.joinpath(
        protein_name, conformation, forcefield, top_name
    )
    traj = md.load(traj_path, top=top_path)
    if t0 is not None:
        # Optionally discard frames before t0 (in ns) to avoid non-equilibrated data.
        traj = traj[traj.time > (t0 * 1e3)]  # Convert t0 from ns to ps

    return traj


def load_spin_data(protein_name: str, type: str) -> pd.DataFrame:
    """Loads experimental spin relaxation data from a CSV file."""
    filename = f"{type}.csv"
    file_path = Path(__file__).parent.joinpath(protein_name, filename)
    df = pd.read_csv(file_path, header=None).rename(columns={0: "res", 1: "value"})
    # Due to uncertain placement of points in digitized plots, need to correct the residue column by rounding to the nearest integer
    df["res"] = df["res"].round().astype(int)
    return df


def find_nh_bonds(traj: md.Trajectory) -> tuple[np.ndarray, np.ndarray]:
    """Return bonded backbone N-H atom indices and their residue numbers."""
    bonds = []
    residues = []
    for atom_a, atom_b in traj.topology.bonds:
        if {atom_a.name, atom_b.name} != {"N", "H"}:
            continue
        # Skip free amino terminus
        if atom_a.residue.resSeq == 1 or atom_b.residue.resSeq == 1:
            continue

        nitrogen = atom_a if atom_a.name == "N" else atom_b
        hydrogen = atom_b if atom_a.name == "N" else atom_a
        bonds.append((nitrogen.index, hydrogen.index))
        residues.append(nitrogen.residue.resSeq)
    if not bonds:
        raise ValueError("The topology contains no bonded N-H pairs.")
    return np.asarray(bonds, dtype=int), np.asarray(residues, dtype=int)


def nh_correlation(
    traj: md.Trajectory,
    nh_bonds: np.ndarray | None = None,
    max_lag: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Calculate residue-resolved second-order N-H rotational correlations.

    MDTraj stores coordinates in nm. Returned lags are in seconds and the
    correlation array has shape (n_lags, n_bonds).
    """
    if traj.n_frames < 2:
        raise ValueError("At least two trajectory frames are required.")
    if nh_bonds is None:
        nh_bonds, _ = find_nh_bonds(traj)
    if max_lag is None:
        max_lag = traj.n_frames // 2
    max_lag = min(max_lag, traj.n_frames - 1)
    if max_lag < 1:
        raise ValueError("max_lag must be at least one frame.")

    vectors = traj.xyz[:, nh_bonds[:, 1], :] - traj.xyz[:, nh_bonds[:, 0], :]
    vectors /= np.linalg.norm(vectors, axis=2, keepdims=True)
    n_frames = vectors.shape[0]
    nfft = 2 * n_frames
    correlations = np.zeros((max_lag, len(nh_bonds)), dtype=float)

    # P2(u.v) = 1.5 * sum_ij(u_i u_j)(v_i v_j) - 0.5.
    for i in range(3):
        for j in range(3):
            component = vectors[:, :, i] * vectors[:, :, j]
            spectrum = np.fft.rfft(component, n=nfft, axis=0)
            autocorrelation = np.fft.irfft(spectrum * spectrum.conj(), n=nfft, axis=0)
            autocorrelation = autocorrelation[:max_lag]
            autocorrelation /= np.arange(n_frames, n_frames - max_lag, -1)[:, None]
            correlations += 1.5 * autocorrelation.real
    correlations -= 0.5
    lags_ps = np.arange(max_lag) * np.mean(np.diff(traj.time))
    # Convert to SI units (seconds) for the output.
    return lags_ps * 1.0e-12, correlations


def fit_correlation(
    lags: np.ndarray,
    correlations: np.ndarray,
    tau_min: float = 1.0e-12,
    tau_max: float = 100.0e-9,
    n_timescales: int = 100,
    fit_max_time: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit each correlation to non-negative exponential weights."""
    if fit_max_time is None:
        fit_max_time = 10.0 * tau_max
    mask = lags <= fit_max_time
    fit_lags = lags[mask]
    if fit_lags.size < 2:
        raise ValueError("The trajectory does not contain enough lag points to fit.")
    taus = np.geomspace(tau_min, tau_max, n_timescales)
    basis = np.exp(-fit_lags[:, None] / taus[None, :])
    weights = np.empty((correlations.shape[1], n_timescales), dtype=float)
    for bond_index in range(correlations.shape[1]):
        weights[bond_index], _ = optimize.nnls(basis, correlations[mask, bond_index])
    return taus, weights


def spectral_density(
    omega: np.ndarray | float, taus: np.ndarray, weights: np.ndarray
) -> np.ndarray:
    """Evaluate J(omega) for one or more fitted correlation functions."""
    omega = np.atleast_1d(np.asarray(omega, dtype=float))
    return 2.0 * np.sum(
        weights[..., None]
        * taus[None, :, None]
        / (1.0 + (omega[None, None, :] * taus[None, :, None]) ** 2),
        axis=1,
    )


def calculate_relaxation(
    taus: np.ndarray,
    weights: np.ndarray,
    device_frequency_mhz: float = 600.0,
) -> pd.DataFrame:
    """Calculate R1, R2, T1, T2, and heteronuclear NOE per N-H bond."""
    # Convert frequency from MHz to magnetic field strength (B0) in units of Tesla
    b0 = 2.0 * np.pi * device_frequency_mhz * 1e6 / GAMMA_H
    omega_h = GAMMA_H * b0
    omega_n = GAMMA_N * b0
    dipolar = MU0 * HBAR * GAMMA_H * GAMMA_N / (4.0 * np.pi * NH_DISTANCE**3)
    frequencies = np.array(
        [0.0, abs(omega_h - omega_n), omega_n, omega_h, omega_h + omega_n]
    )
    j_zero, j_diff, j_n, j_h, j_sum = spectral_density(frequencies, taus, weights).T

    r1 = (dipolar**2 / 20.0) * (j_diff + 3.0 * j_n + 6.0 * j_sum)
    r1 += (CSA_N * omega_n) ** 2 / 15.0 * j_n
    r2 = (dipolar**2 / 40.0) * (
        4.0 * j_zero + j_diff + 3.0 * j_n + 6.0 * j_h + 6.0 * j_sum
    )
    r2 += (CSA_N * omega_n) ** 2 / 90.0 * (4.0 * j_zero + 3.0 * j_n)
    t1 = 1.0 / r1
    t2 = 1.0 / r2
    hetnoe = 1.0 + (dipolar**2 / 20.0) * (6.0 * j_sum - j_diff) * t1 * GAMMA_H / GAMMA_N
    return pd.DataFrame({"R1": r1, "R2": r2, "T1": t1, "T2": t2, "hetNOE": hetnoe})


def save_relaxation_results(
    output_dir: str | Path,
    residues: np.ndarray,
    predictions: pd.DataFrame,
    lags: np.ndarray,
    correlations: np.ndarray,
    taus: np.ndarray,
    weights: np.ndarray,
) -> None:
    """Save predictions as CSV and reusable correlation/fitting data as NPZ."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    predictions = predictions.copy()
    predictions.insert(0, "res", residues)
    predictions.to_csv(output_dir / "spin_relaxation.csv", index=False)
    np.savez_compressed(
        output_dir / "spin_relaxation_intermediate.npz",
        residues=residues,
        lags=lags,
        correlations=correlations,
        taus=taus,
        weights=weights,
    )


def calculate_from_trajectory(
    traj: md.Trajectory,
    output_dir: str | Path | None = None,
    device_frequency_mhz: float = 600.0,
) -> pd.DataFrame:
    """Run the full trajectory-to-relaxation pipeline."""
    bonds, residues = find_nh_bonds(traj)
    lags, correlations = nh_correlation(traj, bonds)
    taus, weights = fit_correlation(lags, correlations)
    predictions = calculate_relaxation(taus, weights, device_frequency_mhz)
    predictions.insert(0, "res", residues)
    if output_dir is not None:
        save_relaxation_results(
            output_dir,
            residues,
            predictions.drop(columns="res"),
            lags,
            correlations,
            taus,
            weights,
        )
    return predictions


devices_mhz = {
    "asyn": 600.0,
    "p53n": 850.0,
    "prta": 600.0,
}


# Device frequencies according to dataset:
# - asyn: 600 MHz "All NMR data were collected at 10 °C using a Varian Unity INOVA 600-MHz spectrometer" https://doi.org/10.1074/jbc.M106777200
# - p53-NTD: 850 MHz "Bruker AVANCE III 850 MHz" https://doi.org/10.13018/BMR52703
# - prothymosin alpha (2013 & 2007 data): "All NMR experiments (except the amide HX experiments)
# were performed on a Varian Inova 600-MHz spectrometer" http://dx.doi.org/10.1016/j.jmb.2013.01.005
# "NMR experiments described below, unless stated otherwise, were performed at 22 °C on a Varian Inova 600 MHz spectrometer equipped with an xyz-gradient triple resonance probe" https://doi.org/10.1021/bi7014822
def main():
    t0_data = load_t0()
    protein = "p53n"
    conf = "compact"
    ff = "c36m"

    frequency = devices_mhz[protein]
    t0 = t0_data[protein][conf][ff]
    traj = load_traj(protein, conf, ff, t0=t0)
    print(traj.n_frames)

    predictions = calculate_from_trajectory(
        traj,
        device_frequency_mhz=frequency,
        output_dir=Path("results") / protein / conf / ff,
    )
    print(predictions.head())


if __name__ == "__main__":
    main()
