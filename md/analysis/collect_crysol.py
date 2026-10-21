from pathlib import Path

import numpy as np


def load_int_file(path: Path, usecols: tuple = (1,)) -> np.ndarray:
    """
    Load a .int file and return the data as a numpy array.
    """
    # Skip the header line (first row)
    # Use only the specified columns; we only need to fetch the first column once, because it is repeeated across all files
    return np.loadtxt(path, skiprows=1, usecols=usecols)


def main():
    root_dir = Path.cwd()
    protein_dirs = [d for d in root_dir.iterdir() if d.is_dir()]
    for protein_dir in protein_dirs:
        print(f"[INFO] Processing {protein_dir.name}...")
        confs = [d for d in protein_dir.iterdir() if d.is_dir()]
        for conf in confs:
            print(f"[INFO] Processing {conf.name}...")
            ffs = [d for d in conf.iterdir() if d.is_dir()]
            for ff in ffs:
                print(f"[INFO] Processing {ff.name}...")
                output_file_name = f"{protein_dir.name}_{conf.name}_{ff.name}.int"
                int_files = sorted(Path(ff).glob("*.int"))
                if not int_files:
                    print(f"[WARNING] No .int files found in {ff.name}. Skipping.")
                    continue
                q_array = load_int_file(int_files[0], usecols=(0,))
                i_array = []
                for int_file in int_files:
                    i_array.append(load_int_file(int_file, usecols=(1,)))
                # array shape N_frames, N_q_points
                i_matrix = np.vstack(i_array)

                np.savez_compressed(output_file_name, q=q_array, i=i_matrix)
                print("[SUCCESS] Saved data to: ", output_file_name)


if __name__ == "__main__":
    main()
