import pandas as pd
import subprocess
import os

df = pd.read_csv("ligands_clean_admet_pass.csv")

os.makedirs("pdbqt", exist_ok=True)

for _, row in df.iterrows():

    cid = str(row["Compound_ID"])
    smiles = str(row["SMILES"])

    command = [
        "obabel",
        f"-:{smiles}",
        "--gen3d",
        "-O",
        f"pdbqt/{cid}.pdbqt"
    ]

    subprocess.run(command)

print("Finished!")
