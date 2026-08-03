#!/usr/bin/env python3
"""
filter_admet.py

Filters a library of ligands by computed ADMET-relevant molecular properties.

This script uses RDKit to compute physicochemical descriptors that are
commonly used as ADMET (Absorption, Distribution, Metabolism, Excretion,
Toxicity) proxies. It is NOT a substitute for a real ADMET prediction tool
(SwissADME, pkCSM, ADMETlab, etc.) -- those use trained ML models for things
like CYP450 inhibition, hERG liability, or hepatotoxicity that cannot be
derived from structure alone with simple descriptors. What this script DOES
give you is a fast, reproducible, fully offline first-pass filter on the
physicochemical properties that drive most ADME behavior:

    - Molecular weight (MW)              -- absorption / permeability
    - LogP                                -- lipophilicity, membrane permeation
    - TPSA (topological polar surf. area)-- absorption, blood-brain barrier
    - H-bond donors / acceptors           -- absorption, permeability
    - Rotatable bonds                     -- oral bioavailability (flexibility)
    - Number of aromatic rings            -- metabolic stability / promiscuity risk
    - Molar refractivity                  -- general "drug-likeness" (Lipinski MR rule)

Default thresholds implement a combination of:
    - Lipinski's Rule of Five (MW<=500, LogP<=5, HBD<=5, HBA<=10)
    - Veber's rules for oral bioavailability (TPSA<=140, RotBonds<=10)

You can override every threshold from the command line.

Usage:
    python filter_admet.py input.csv --smiles-col SMILES
    python filter_admet.py input.csv --smiles-col SMILES --max-logp 4.5 --max-tpsa 120
"""

import argparse
import sys

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from rdkit.RDLogger import DisableLog

DisableLog("rdApp.*")


def compute_admet_properties(smiles: str) -> dict | None:
    """Compute ADMET-proxy descriptors for a single SMILES string.

    Returns None if RDKit cannot parse the SMILES.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return {
        "MW": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "RotBonds": Descriptors.NumRotatableBonds(mol),
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "MolarRefractivity": Crippen.MolMR(mol),
        "HeavyAtoms": mol.GetNumHeavyAtoms(),
    }


def passes_filters(props: dict, args) -> tuple[bool, list[str]]:
    """Check a property dict against all thresholds.

    Returns (passed, list_of_violated_rule_names).
    """
    violations = []

    if props["MW"] > args.max_mw:
        violations.append(f"MW>{args.max_mw}")
    if props["LogP"] > args.max_logp:
        violations.append(f"LogP>{args.max_logp}")
    if props["LogP"] < args.min_logp:
        violations.append(f"LogP<{args.min_logp}")
    if props["TPSA"] > args.max_tpsa:
        violations.append(f"TPSA>{args.max_tpsa}")
    if props["HBD"] > args.max_hbd:
        violations.append(f"HBD>{args.max_hbd}")
    if props["HBA"] > args.max_hba:
        violations.append(f"HBA>{args.max_hba}")
    if props["RotBonds"] > args.max_rotbonds:
        violations.append(f"RotBonds>{args.max_rotbonds}")
    if props["AromaticRings"] > args.max_aromatic_rings:
        violations.append(f"AromaticRings>{args.max_aromatic_rings}")

    passed = len(violations) <= args.max_violations
    return passed, violations


def main():
    parser = argparse.ArgumentParser(
        description="Filter ligands by RDKit-computed ADMET-relevant properties "
                    "(Lipinski + Veber rules by default)."
    )
    parser.add_argument("input_csv", help="Path to input CSV file containing ligands.")
    parser.add_argument(
        "--smiles-col", default="SMILES",
        help="Name of the column containing SMILES strings (default: SMILES)."
    )

    # Thresholds -- defaults combine Lipinski (Ro5) + Veber rules
    parser.add_argument("--max-mw", type=float, default=500.0, help="Max molecular weight (default: 500)")
    parser.add_argument("--max-logp", type=float, default=5.0, help="Max LogP (default: 5)")
    parser.add_argument("--min-logp", type=float, default=-2.0, help="Min LogP (default: -2; too hydrophilic also hurts permeability)")
    parser.add_argument("--max-tpsa", type=float, default=140.0, help="Max TPSA in A^2 (default: 140, Veber rule)")
    parser.add_argument("--max-hbd", type=int, default=5, help="Max H-bond donors (default: 5)")
    parser.add_argument("--max-hba", type=int, default=10, help="Max H-bond acceptors (default: 10)")
    parser.add_argument("--max-rotbonds", type=int, default=10, help="Max rotatable bonds (default: 10, Veber rule)")
    parser.add_argument("--max-aromatic-rings", type=int, default=4, help="Max aromatic rings (default: 4; >4 associated with poor solubility/promiscuity)")
    parser.add_argument(
        "--max-violations", type=int, default=1,
        help="Maximum number of rule violations allowed to pass (default: 1)."
    )

    parser.add_argument("--output", default=None, help="Output CSV for passing compounds.")
    parser.add_argument("--rejected-output", default=None, help="Output CSV for rejected compounds.")
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

    prop_cols = {
        "MW": [], "LogP": [], "TPSA": [], "HBD": [], "HBA": [],
        "RotBonds": [], "AromaticRings": [], "MolarRefractivity": [], "HeavyAtoms": [],
    }
    valid_list = []
    pass_list = []
    violation_list = []

    n_invalid = 0
    for smi in df[args.smiles_col]:
        smi_str = str(smi) if pd.notna(smi) else ""
        props = compute_admet_properties(smi_str)

        if props is None:
            n_invalid += 1
            for k in prop_cols:
                prop_cols[k].append(None)
            valid_list.append(False)
            pass_list.append(False)
            violation_list.append("INVALID_SMILES")
            continue

        for k in prop_cols:
            prop_cols[k].append(props[k])

        passed, viol = passes_filters(props, args)
        valid_list.append(True)
        pass_list.append(passed)
        violation_list.append(";".join(viol) if viol else "")

    for k, v in prop_cols.items():
        df[k] = v
    df["Valid_SMILES"] = valid_list
    df["ADMET_Pass"] = pass_list
    df["Violations"] = violation_list

    if n_invalid > 0:
        print(
            f"WARNING: {n_invalid} of {len(df)} rows had SMILES that RDKit could not parse "
            f"and were excluded from the pass set. Check your source data.",
            file=sys.stderr,
        )

    passed_df = df[df["ADMET_Pass"]].copy()
    rejected_df = df[~df["ADMET_Pass"]].copy()

    out_path = args.output or args.input_csv.rsplit(".", 1)[0] + "_admet_pass.csv"
    rej_path = args.rejected_output or args.input_csv.rsplit(".", 1)[0] + "_admet_rejected.csv"

    passed_df.to_csv(out_path, index=False)
    rejected_df.to_csv(rej_path, index=False)

    print(f"\nTotal compounds:   {len(df)}")
    print(f"Invalid SMILES:    {n_invalid}")
    print(f"Passed:            {len(passed_df)}")
    print(f"Rejected:          {len(rejected_df)}")
    print(f"\nPassed compounds written to:   {out_path}")
    print(f"Rejected compounds written to: {rej_path}")


if __name__ == "__main__":
    main()
