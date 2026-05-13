# Project Titlte: Molecular Docking and MD Simulation-Driven Discovery of Zn2+ Targeting MMP-1 Inhibitors for Next-Generation Bioactive Sunscreens 

## Project Summary
Chronic UV exposure triggers the overproduction of Matrix Metalloproteinase-1 (MMP-1), a collagenase that initiates the degradation of the dermal extracellular matrix. Our project studies the ability for zinc-binding ligands to inhibit MMP-1. The computational workflow uses molecular docking simulations to screen for the top 3 candidate ligands from a known library, followed by molecular dynamics to simulate ligand stability in the binding pocket over time. The main output is a ranking of the ideal inhibitor molecules that could potentially be implemented to support next-generation sunscreens.

## Package Contents

- `Molecular Docking` folder
  - `ligands`
  - `ligands_raw`
  - `results`
  - `table_results`
  - `1cgl_no_lig.pdb`
  - `1CGL.PDB`
  - `config.txt`
  - `environment.yml` (for docking)
  - `flex_receptor_flex.pdbqt`
  - `flex_receptor_rigid.pdbqt`
  - `flex_receptor.box.txt`
  - `flex_receptor.box.txt`
  - `Lab5-Docking.md`
  - `MMP1_receptor.pdb`
  - `Molecular Docking.ipynb`
  - `Molecular Docking(Dan ver).ipynb`
- `Molecular Dynamics` folder
  - `assets` (contains 1CGL; ligands 0 (control), 2, 8, 9)
  - `md_input` (json files for running md simulations on ligands)
  - `md_results` (md data on tested ligands)
  - `references` (lab notebooks that were refered to)
  - `environment.yml` (for molecular dynamics)
  - `md.ipynb` (notebook to generate md simulations and RMSD vs Time on data)
    
## Environment Setup

For Molecular Docking (navigate into Molecular Docking folder):
Create the environment from YAML:
- `conda env create -f environment.yml`
- `conda activate vina`

For Molecular Dynamics (MD) (navigate into Molecular Dynamics folder):
Create the environment from YAML:
- `conda env create -f environment.yml`
- `conda activate md`

## How To Run

(docking notebook)

For Molecular Dynamics (open md.ipynb):
Step 1: Run the PDBFixer cell to fix assets/1cgl_no_lig.pdb (adds missing hydrogens). Outputs assets/1cgl_fixed.pdb.
Step 2: Run subsequent cells to perform MD simulations (define vacuum ChemicalSystem for ligands, set simulation settings [1 ns production run], create Protocol & NonTransformation, and write JSON files to run in termonal)
- NOTE: system is performed in vacuum to improve runtime. In creating the settings, the settings for running NVT and NPT are purposely set to be really low to reduce runtime (our MD study only concerned about RMSD).
Step 3. Run RMSD cell to generate a plot of Backbone RMSD vs. Time, along with average and max RMSD values for each ligand

## Expected Outputs

For Molecular Dynamics:
- MAIN OUTPUT `rmsd.png` Backbone RMSD vs. Time plot
- `md_results` (each ligand folder contains: equil_npt.pdb, equil_nvt.pdb, minimized.pdb, simulation.log)

## Runtime Notes And Limitations

For MD, each ligand takes ~30-60 minutes on CPU (for a 1ns simulation). For our purposes, a 1 ns runtime is mostly sufficient (RMSD results plateus and mostly stabilizes). Though a GPU is not required to run, it would greatly speed up the runtime of the simulations. A vacuum environment is created to reduce runtime, and so RMSD values may not accurately reflect how the ligand wound perform in physiological solvent conditions. 
