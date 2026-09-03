import mdtraj as md
import numpy as np
import pandas as pd


def compute(
    dir: str,
    comp_rg: bool = True,
    comp_cs: bool = False,
    comp_ss: bool = True,
    xtc_file="md_protein.xtc",
    top_file="protein.gro",
):
    """
    Compute radius of gyration and/or chemical shifts from the trajectory loaded from the specified directory.
    """
    # Load trajectory
    traj = md.load_xtc(f"{dir}/{xtc_file}", top=f"{dir}/{top_file}")
    # Select protein atoms
    protein = traj.topology.select("protein")
    protein_traj = traj.atom_slice(protein)
    if comp_rg:
        print("Getting radius of gyration...")
        # Get masses of atoms in protein
        massess = np.array([atom.element.mass for atom in protein_traj.topology.atoms])
        rg = md.compute_rg(protein_traj, masses=massess) * 10  # convert to angstrom
        time = protein_traj.time / 1000  # convert to ns
        # Save as CSV
        rg_df = pd.DataFrame({"time": time, "rg": rg})
        rg_df.to_csv(f"{dir}/rg.csv", index=False)
    if comp_cs:
        n_frames = protein_traj.n_frames
        chunk_size = n_frames // 3
        for i, start in enumerate(range(0, n_frames, chunk_size)):
            end = min(start + chunk_size, n_frames)
            print(
                f"Processing chunk {i + 1}/{(n_frames + chunk_size - 1) // chunk_size} from frame {start} to {end}"
            )
            # Slice the trajectory for the current chunk
            chunk_traj = protein_traj[start:end]
            # Compute chemical shifts for the current chunk
            shifts = md.chemical_shifts_spartaplus(chunk_traj, rename_HN=True)
            # Save the shifts to a pickle file
            shifts.columns = chunk_traj.time / 1000
            shifts.to_pickle(f"{dir}/shifts_{i}.pkl")
    if comp_ss:
        dssp = md.compute_dssp(protein_traj, simplified=True)
        np.save(f"{dir}/dssp.npy", dssp)
        print("Saved DSSP data to dssp.npy")


# Example usage:
# Note - The chemical shifts computation requires the SPARTA+ software to be installed and properly configured.
def main():
    protein_name = "p53n"
    conformation = "extended"
    forcefield = "a99disp"
    compute(
        dir=f"{protein_name}/{conformation}/{forcefield}",
        comp_rg=True,
        comp_cs=True,
        comp_ss=True,
        xtc_file="md_protein.xtc",
        top_file="protein.gro",
    )


if __name__ == "__main__":
    main()
