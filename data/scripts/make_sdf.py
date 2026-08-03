import pandas as pd
import subprocess

df = pd.read_csv("ligands_clean_admet_pass.csv")

selected = [
    "MMV1782354",
    "MMV1581559",
    "MMV1580850"
]

for _, row in df[df["Compound_ID"].isin(selected)].iterrows():
    ligand = row["Compound_ID"]
    smiles = row["SMILES"]

    cmd = [
        "obabel",
        f"-:{smiles}",
        "-O",
        f"{ligand}.sdf",
        "--gen3d"
    ]

    subprocess.run(cmd)
    print(f"Generated {ligand}.sdf")
