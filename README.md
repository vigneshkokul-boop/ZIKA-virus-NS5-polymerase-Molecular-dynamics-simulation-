# Zika Virus NS5 Polymerase Molecular Dynamics Simulation

> **An end-to-end computational drug discovery workflow for identifying potential inhibitors of the Zika virus NS5 RNA-dependent RNA polymerase (RdRp) using molecular docking and molecular dynamics simulations.**

---

## Project Overview

This repository contains the complete computational workflow used to investigate potential antiviral compounds targeting the **Zika virus NS5 RNA-dependent RNA polymerase (RdRp)**. The study combines **protein preparation, ligand library processing, virtual screening, molecular docking, molecular dynamics simulations, and trajectory analysis** to evaluate the stability and binding behavior of candidate inhibitors.

The project demonstrates a complete **structure-based drug discovery (SBDD)** pipeline using open-source computational biology and cheminformatics tools.

---

## Project Status

🟢 **Status:** Ongoing

This project is being carried out as part of an **ongoing Bioinformatics Internship** at the **Rajiv Gandhi Centre for Biotechnology (RGCB), Thiruvananthapuram, India**, under the guidance of **Dr. Kathiresan**.

The repository will continue to be updated as additional analyses, scripts, and simulation results become available.

---

# Why is this project important?

Zika virus is a mosquito-borne flavivirus capable of causing severe neurological complications, including congenital Zika syndrome and Guillain–Barré syndrome. Despite extensive research, there are currently no approved antiviral drugs that specifically target the viral RNA-dependent RNA polymerase.

Computational drug discovery techniques such as **virtual screening**, **molecular docking**, and **molecular dynamics (MD) simulations** enable the rapid identification and evaluation of potential inhibitors before experimental validation, significantly reducing both time and cost in the drug discovery process.

---

# Objectives

- Retrieve and prepare the Zika NS5 polymerase structure.
- Prepare the MMV compound library for virtual screening.
- Filter compounds using Lipinski's Rule of Five and ADME properties.
- Perform molecular docking using AutoDock Vina.
- Identify promising lead compounds.
- Validate protein–ligand complexes using molecular dynamics simulations.
- Analyze structural stability and conformational dynamics.
- Identify potential antiviral candidates for future experimental studies.

---

# Methodology

## 1. Protein Preparation

- Crystal structure retrieved from the **Protein Data Bank (PDB ID: 5WZ3)**.
- Missing loop regions modeled and refined.
- Protein prepared using **CHARMM-GUI PDB Reader & Manipulator**.
- Final processed protein generated as:

```
step1_pdbreader.pdb
```

---

## 2. Ligand Library Preparation

The ligand library was obtained from the **Medicines for Malaria Venture (MMV)** compound database.

Workflow:

```
MMV Excel Dataset
        │
        ▼
LibreOffice
(Excel → CSV)
        │
        ▼
CSV → SMILES
(Python + Open Babel)
        │
        ▼
Lipinski Filtering
        │
        ▼
ADME Filtering
        │
        ▼
RDKit
(3D Structure Generation)
        │
        ▼
SDF Ligands
```

---

## 3. Molecular Docking

Protein receptor preparation was performed using AutoDock Tools.

Docking simulations were carried out using **AutoDock Vina** to predict:

- Binding affinity
- Binding pose
- Protein–ligand interactions

Top-ranked compounds were selected for molecular dynamics simulations.

---

## 4. Molecular Dynamics Simulation

Selected protein–ligand complexes were subjected to molecular dynamics simulations using **GROMACS**.

Trajectory analyses included:

- Backbone RMSD
- RMSF
- Radius of Gyration (Rg)
- Hydrogen Bond Analysis
- Principal Component Analysis (PCA)
- MM/PBSA Binding Free Energy

---

# Python Scripts

This repository includes custom Python scripts for automating ligand preparation and virtual screening.

Included scripts:

- CSV to SMILES conversion
- Lipinski Rule of Five filtering
- ADME filtering
- RDKit 3D ligand generation
- Receptor preparation

---

# Repository Structure

```
ZIKA-virus-NS5-polymerase-Molecular-dynamics-simulation/

├── data/
│   ├── protein/
│   ├── ligands/
│   └── README.md
│
├── scripts/
│   ├── csv_to_smiles.py
│   ├── lipinski_filter.py
│   ├── adme_filter.py
│   ├── rdkit_3d_conversion.py
│   ├── prepare_receptor.py
│   └── README.md
│
├── results/
│   ├── RMSD
│   ├── RMSF
│   ├── PCA
│   ├── Docking Results
│   └── README.md
│
├── Thesis/
│
├── images/
│
├── README.md
└── LICENSE
```

---

# Software and Tools

- GROMACS
- CHARMM-GUI
- AutoDock Vina
- Open Babel
- RDKit
- Python
- NumPy
- Pandas
- Matplotlib
- MDAnalysis
- LibreOffice

---

# Key Outcomes

- Successfully prepared the Zika NS5 polymerase structure for molecular docking and molecular dynamics simulations.
- Processed the MMV compound library into docking-ready ligand structures.
- Filtered compounds using Lipinski and ADME criteria.
- Performed molecular docking to identify promising lead compounds.
- Conducted molecular dynamics simulations to evaluate protein–ligand stability.
- Generated publication-quality trajectory analyses including RMSD, RMSF, PCA, and binding free energy calculations.

---

# Reproducibility

To reproduce this workflow:

1. Download the protein structure (PDB ID: **5WZ3**).
2. Prepare the protein using CHARMM-GUI.
3. Download the MMV compound library.
4. Convert the ligand dataset to SMILES.
5. Apply Lipinski and ADME filtering.
6. Generate 3D ligand structures using RDKit.
7. Perform molecular docking using AutoDock Vina.
8. Select top compounds for molecular dynamics simulations in GROMACS.
9. Analyze trajectories using the provided Python scripts.

---

# Future Work

- Extend simulations to additional lead compounds.
- Perform replicate MD simulations.
- Compare multiple antiviral compounds.
- Develop machine learning models for binding affinity prediction.
- Validate computational findings through experimental collaboration.

---

# Acknowledgements

This work is being carried out as part of an **ongoing Bioinformatics Internship** at the **Rajiv Gandhi Centre for Biotechnology (RGCB), Thiruvananthapuram, India**.

I sincerely thank **Dr. Kathiresan Natarajan** for his guidance, mentorship, and support throughout this project.

---

# Author

Kokulraj A B 

**M.Sc. Food tech**

**Research Interests**

- Structural Bioinformatics
- Molecular Dynamics Simulation
- Structure-Based Drug Discovery
- Computational Biology
- Machine Learning for Drug Discovery

---

## License

This project is released under the **MIT License**.
