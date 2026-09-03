from io import load_dssp, load_rg, process_dssp

import numpy as np

TAU_MUL = 3.0


def autocorrelation(array: np.ndarray, max_lag: int | None = None) -> np.ndarray:
    x = array - np.mean(array)
    n = len(x)
    if max_lag is None:
        max_lag = n // 2
    # fft-based correlation
    f = np.fft.fft(x, n * 2)
    acf = np.fft.ifft(f * np.conjugate(f))[:n].real
    # normalization
    acf /= acf[0]
    return acf[:max_lag]


def integrated_tau_running(lags: np.ndarray, acf: np.ndarray) -> float:
    dt = np.mean(np.diff(lags))
    tau_run = np.cumsum(acf) * dt
    return max(tau_run)


def get_tau(
    time: np.ndarray,
    x: np.ndarray,
    max_lag: int | None = None,
):
    # simulation timestep
    dt = np.mean(np.diff(time))
    acf = autocorrelation(x, max_lag=max_lag)
    lags = np.arange(len(acf)) * dt

    tau = integrated_tau_running(lags, acf)
    return lags, acf, tau


def exp_decay(x: np.ndarray, tau: float):
    """
    Simple exponential decay function representing an approximation of `C(t)` with `A=1.0`.
    """
    # A = 1.0 forces the function to pass through y = 1.0 when x = 0
    return 1.0 * np.exp(-x / tau)


def main():
    protein_name = "p53n"
    conformation = "extended"
    forcefield = "a99disp"
    rg_df = load_rg(f"{protein_name}/{conformation}/{forcefield}/rg.csv")
    dssp_df = load_dssp(
        f"{protein_name}/{conformation}/{forcefield}/dssp.npy",
        time_array=np.array(rg_df["time"].values),
    )
    helix_df = process_dssp(dssp_df, code="H", timeseries=True)
    beta_df = process_dssp(dssp_df, code="E", timeseries=True)

    rg_lags, rg_acf, rg_tau = get_tau(rg_df["time"], rg_df["rg"])
    helix_lags, helix_acf, helix_tau = get_tau(helix_df["time"], helix_df["fraction"])
    beta_lags, beta_acf, beta_tau = get_tau(beta_df["time"], beta_df["fraction"])

    data = {
        "rg": {"tau": rg_tau, "lags": rg_lags, "acf": rg_acf},
        "helix": {"tau": helix_tau, "lags": helix_lags, "acf": helix_acf},
        "beta": {"tau": beta_tau, "lags": beta_lags, "acf": beta_acf},
    }

    max_tau = max(data[key]["tau"] for key in data)
    print(max_tau)
    t0 = max_tau * 3.0
    print(t0)  # ns
