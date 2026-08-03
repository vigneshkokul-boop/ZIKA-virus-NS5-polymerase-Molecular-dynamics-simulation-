from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import os

# Read CSV
df = pd.read_csv("compounds_admet_pass.csv", encoding="latin-1")

# Create output folder
os.makedirs("sdf_files", exist_ok=True)

success = 0
failed = 0

for _, row in df.iterrows():

    smiles = str(row["SMILES"]).strip()
    mmv = str(row["MMV ID"]).strip()

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        failed += 1
        continue

    mol = Chem.AddHs(mol)

    status = AllChem.EmbedMolecule(
        mol,
        AllChem.ETKDGv3()
    )

    if status != 0:
        failed += 1
        continue

    AllChem.UFFOptimizeMolecule(mol)

    mol.SetProp("_Name", mmv)

    writer = Chem.SDWriter(f"sdf_files/{mmv}.sdf")
    writer.write(mol)
    writer.close()

    success += 1

print(f"Generated SDF files : {success}")
print(f"Failed              : {failed}")
