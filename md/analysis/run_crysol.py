import subprocess
from contextlib import chdir
from pathlib import Path

from tqdm import tqdm

LM = 50
FB = 18


def run_crysol_cmd(file_name: str) -> None:
    cmd = f"crysol --explicit-hydrogens --fb {FB} --lm {LM} {file_name}"
    res = subprocess.run(cmd, shell=True, check=True, encoding="utf-8")
    if res.returncode == 0:
        print(f"[INFO] Crysol completed successfully for {file_name}")
    else:
        print(f"[ERROR] Crysol failed for {file_name}. Message: {res.stderr}")


def main():
    # Run this from within the pdbs directory
    root_dir = Path.cwd()
    pdb_files = root_dir.rglob("*.pdb")

    for file_path in tqdm(pdb_files, desc="Processing PDB files"):
        # Change working directory to the directory containing the PDB file, so that outputs are saved to this directory
        with chdir(file_path.parent):
            run_crysol_cmd(file_path.name)


if __name__ == "__main__":
    main()
