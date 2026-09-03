import numpy as np

BLOCK_MUL = 5.0
# Block size is BLOCK_MUL * tau_c


def get_x_prod(t0: dict, data: dict) -> dict:
    """
    data: dictionary with conformations as keys and time, variable of interest dataframe as values
    Filters "X" in each trajectory based on the autocorrelation time
    """
    x_prod = {}
    for conf, item in data.items():
        # Equilibration time
        if item.__class__ == dict:
            # Chemical shifts
            x_prod[conf] = {}
            for shift, df in item.items():
                current_t0 = t0[conf]
                mask = df.time >= current_t0
                df = df[mask].reset_index(drop=True)
                x_prod[conf][shift] = df
        else:
            current_t0 = t0[conf]
            mask = item.time >= current_t0
            df = item[mask].reset_index(drop=True)
            x_prod[conf] = df
    return x_prod


def split_blocks(x: np.array, block_size: int, squared: bool = False) -> list:
    """
    Splits data into blocks based on the autocorrelation time.  Squared = True for Rg
    """
    num_frames = len(x)
    num_blocks = num_frames // block_size
    blocks = []
    for i in range(num_blocks + 1):
        start = i * block_size
        end = start + block_size
        # Discard last block if number of elements is less than block size
        if end > num_frames:
            break
        block = x[start:end] ** (2 if squared else 1)
        blocks.append(block)
    return blocks


def block_averaging(
    blocks: list, weights: list | None = None, squared: bool = False
) -> tuple:
    """
    Given a list of blocks, computes the mean and standard error of the mean of the block means.  Returns (mean, se).
    Formula for unbiased estimator of the standard error of the weighted mean taken from https://seismo.berkeley.edu/~kirchner/Toolkits/Toolkit_12.pdf "Case 1" where we want to give more weight to a set of points
    """
    blocks_means = np.array([np.mean(block) for block in blocks])
    if weights is None:
        weights = np.array([len(block) for block in blocks])
    # weighted mean = unweighted mean if all weights equal
    block_mean = np.average(blocks_means, weights=weights)
    # neff = length of blocks if all weights equal
    neff = np.sum(weights) ** 2 / np.sum(weights**2)
    block_var = (
        neff
        / (neff - 1)
        * np.sum(weights * (blocks_means - block_mean) ** 2)
        / np.sum(weights)
    )
    block_std = block_var**0.5
    block_se = block_std / np.sqrt(neff)
    if not squared:
        return block_mean, block_se
    # If squared, return sqrt of mean and propagate error
    # block mean => <Rg^2> => square root to get Rg
    block_std = 1.0 / (2.0 * block_mean**0.5) * block_std
    block_se = block_std / np.sqrt(neff)
    return block_mean**0.5, block_se


def block_bootstrapping(
    blocks: list, bin_edges: np.ndarray, nboot: int = 1000
) -> tuple:
    boot_counts = np.zeros((nboot, len(bin_edges) - 1))
    nb = len(blocks)
    for b in range(nboot):
        idxs = np.random.randint(0, nb, size=nb)
        sample = np.concatenate([blocks[i] for i in idxs])
        counts, _ = np.histogram(sample, bins=bin_edges, density=False)
        boot_counts[b] = counts / counts.sum()
    mean_hist = boot_counts.mean(axis=0)
    lower = np.percentile(boot_counts, 2.5, axis=0)
    upper = np.percentile(boot_counts, 97.5, axis=0)
    return mean_hist, lower, upper


def ensemble_average(
    x_prod: dict, taus: dict, col: str, squared: bool = False
) -> tuple:
    """Performs block averaging over concatenated trajectories taking into account the number of independent samples (number of blocks) as weights.
    x_prod: dictionary with conformations as keys and time, variable of interest dataframe as values
    taus: dictionary with conformations as keys corresponding tau_c as values (value specific to the observable of interest)
    col: column name of variable of interest in the dataframes
    squared: if True, computes <X^2> and returns sqrt(<X^2>) with propagated error. If False, computes <X>.
    """
    replica_blocks = []
    replica_weights = []
    conformations = x_prod.keys()

    for conf in conformations:
        # tau_c depends on variable of choice and conformation (replica)
        tau_c = taus[conf]
        block_size = int(np.ceil(BLOCK_MUL * tau_c))
        # get variable of interest as np array
        x = x_prod[conf][col].values
        # split blocks (not overlapping)
        blocks = split_blocks(x, block_size, squared=squared)
        replica_blocks.extend(blocks)
        # weights (same within the same conformation/replica)
        weights = [len(x) / block_size] * len(blocks)
        replica_weights.extend(weights)

    block_mean, block_se = block_averaging(
        replica_blocks, weights=np.array(replica_weights), squared=squared
    )
    return block_mean, block_se, replica_blocks, replica_weights
