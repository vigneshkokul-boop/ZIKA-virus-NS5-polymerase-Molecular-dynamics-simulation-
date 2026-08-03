#!/usr/bin/env python3
"""
filter_lipinski.py

Filters a library of ligands by Lipinski's Rule of Five.

Lipinski's Rule of Five (drug-likeness heuristic for oral bioavailability):
    - Molecular weight       <= 500 Da
    - LogP (octanol-water)   <= 5
    - H-bond donors          <= 5
    - H-bond acceptors       <= 10

A compound is usually considered "Lipinski-compliant" if it violates
AT MOST ONE of these four rules (the standard convention). This script
reports the violation count per compound and lets you set the allowed
number of violations.

Usage:
    python filter_lipinski.py input.csv --smiles-col SMILES --max-violations 1
"""

import argparse
import sys

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
from rdkit.RDLogger import DisableLog

# Silence RDKit's verbose parsing warnings (we handle invalid SMILES ourselves)
DisableLog("rdApp.*")


def compute_lipinski_properties(smiles: str) -> dict | None:
    """Compute the four Lipinski descriptors for a single SMILES string.

    Returns None if the SMILES string could not be parsed by RDKit.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return {
        "MW": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
    }


def count_violations(props: dict) -> int:
    """Count how many of the 4 Lipinski rules are violated."""
    violations = 0
    if props["MW"] > 500:
        violations += 1
    if props["LogP"] > 5:
        violations += 1
    if props["HBD"] > 5:
        violations += 1
    if props["HBA"] > 10:
        violations += 1
    return violations


def main():
    parser = argparse.ArgumentParser(description="Filter ligands by Lipinski's Rule of Five.")
    parser.add_argument("input_csv", help="Path to input CSV file containing ligands.")
    parser.add_argument(
        "--smiles-col", default="SMILES",
        help="Name of the column containing SMILES strings (default: SMILES)."
    )
    parser.add_argument(
        "--max-violations", type=int, default=1,
        help="Maximum number of Lipinski rule violations allowed to pass the filter "
             "(default: 1, the standard convention)."
    )
    parser.add_argument(
        "--output", default=None,
        help="Path for the filtered (passing) compounds CSV. "
             "Default: '<input>_lipinski_pass.csv'"
    )
    parser.add_argument(
        "--rejected-output", default=None,
        help="Path for the rejected compounds CSV, written for audit purposes. "
             "Default: '<input>_lipinski_rejected.csv'"
    )
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input_csv)
    except FileNotFoundError:
        sys.exit(f"Error: file not found: {args.input_csv}")

    if args.smiles_col not in df.columns:
        sys.exit(
            f"Error: column '{args.smiles_col}' not found in {args.input_csv}.\n"
            f"Available columns: {list(df.columns)}"
        )

    mw_list, logp_list, hbd_list, hba_list, viol_list, valid_list = (
        [], [], [], [], [], []
    )

    n_invalid = 0
    for smi in df[args.smiles_col]:
        smi_str = str(smi) if pd.notna(smi) else ""
        props = compute_lipinski_properties(smi_str)

        if props is None:
            n_invalid += 1
            mw_list.append(None)
            logp_list.append(None)
            hbd_list.append(None)
            hba_list.append(None)
            viol_list.append(None)
            valid_list.append(False)
            continue

        mw_list.append(props["MW"])
        logp_list.append(props["LogP"])
        hbd_list.append(props["HBD"])
        hba_list.append(props["HBA"])
        viol_list.append(count_violations(props))
        valid_list.append(True)

    df["MW"] = mw_list
    df["LogP"] = logp_list
    df["HBD"] = hbd_list
    df["HBA"] = hba_list
    df["Lipinski_Violations"] = viol_list
    df["Valid_SMILES"] = valid_list

    if n_invalid > 0:
        print(
            f"WARNING: {n_invalid} of {len(df)} rows had SMILES that RDKit could not parse. "
            f"These are flagged Valid_SMILES=False and excluded from the pass set.\n"
            f"This usually means upstream data corruption (e.g. placeholder/wildcard SMILES) "
            f"-- check your source file before trusting any results below.",
            file=sys.stderr,
        )

    passed = df[(df["Valid_SMILES"]) & (df["Lipinski_Violations"] <= args.max_violations)].copy()
    rejected = df[~((df["Valid_SMILES"]) & (df["Lipinski_Violations"] <= args.max_violations))].copy()

    out_path = args.output or args.input_csv.rsplit(".", 1)[0] + "_lipinski_pass.csv"
    rej_path = args.rejected_output or args.input_csv.rsplit(".", 1)[0] + "_lipinski_rejected.csv"

    passed.to_csv(out_path, index=False)
    rejected.to_csv(rej_path, index=False)

    print(f"\nTotal compounds:        {len(df)}")
    print(f"Invalid SMILES:         {n_invalid}")
    print(f"Passed (<= {args.max_violations} violation(s)):  {len(passed)}")
    print(f"Rejected:               {len(rejected)}")
    print(f"\nPassed compounds written to:   {out_path}")
    print(f"Rejected compounds written to: {rej_path}")


if __name__ == "__main__":
    main()
