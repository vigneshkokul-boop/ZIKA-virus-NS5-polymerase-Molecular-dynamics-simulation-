import csv

with open("compounds.csv", newline="", encoding="latin-1") as infile, \
     open("compounds.smi", "w", encoding="utf-8") as outfile:

    reader = csv.DictReader(infile)

    for row in reader:
        smiles = row["SMILES"].strip()
        mmv = row["MMV ID"].strip()

        if smiles:
            outfile.write(f"{smiles} {mmv}\n")
