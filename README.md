# Project Title: Molecular Docking and MD Simulation-Driven Discovery of Zn2+ Targeting MMP-1 Inhibitors for Next-Generation Bioactive Sunscreens 

## Project Summary
Chronic UV exposure triggers the overproduction of Matrix Metalloproteinase-1 (MMP-1), a collagenase that initiates the degradation of the dermal extracellular matrix. Our project studies the ability for zinc-binding ligands to inhibit MMP-1. The computational workflow uses molecular docking simulations to screen for the top 3 candidate ligands from a known library, followed by molecular dynamics to simulate ligand stability in the binding pocket over time. The main output is a ranking of the ideal inhibitor molecules that could potentially be implemented to support next-generation sunscreens.

## Package Contents

- `Molecular Docking` folder
  - `ligands` (20 candidate ligand pdbqt)
  - `ligands_raw` (sdf files for all ligands)
  - `results` (vina docking output)
  - `table_results` (shows affinity scores)
  - `1cgl_no_lig.pdb` (receptor with ligand removed)
  - `1CGL.PDB` (original receptor from PDB database)
  - `config.txt` (AutoDock Vina config)
  - `environment.yml` (for docking)
  - `flex_receptor_flex.pdbqt` (flexible residues (His 218, 222, 228) of receptor)
  - `flex_receptor_rigid.pdbqt` (rigid portion of prepared receptor)
  - `flex_receptor.box.txt` (defined docking box)
  - `MMP1_receptor.pdb` (receptor after PDBFixer (hydrogens added at pH 7.4))
  - `Molecular Docking.ipynb` (notebook to generate docking simulations)
  - `references` (notebooks that were referred to)
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
- `conda install ipykernel -y python -m ipykernel install --user --name vina --display-name "vina"`

For Molecular Dynamics (MD) (navigate into Molecular Dynamics folder):  
Create the environment from YAML:  
- `conda env create -f environment.yml`
- `conda activate md`

## How To Run

For Molecular Docking (open Molecular Docking.ipynb):  
  
Step 1: Use PDBFixer to clean structure and add missing hydrogens, prepare receptor and ligands. (pH 7.4) (run cells under "For our Project (Just paste all of the rest things in Terminal)")  
Step 2: Prepare all ligands (20)  
Step 3: Run docking for all ligands. This loops over all 20 ligands with AutoDock Vina. (run cells under "IN NOTEBOOK") Poses and log files are saves to `results`.  
Step 4: Rank ligands based ob best binding affinity. The final ligands selected for MD simulations are ligands 2, 8, and 9, along with ligand 0 as a control.  

For Molecular Dynamics (open md.ipynb):  
  
Step 1: Run the PDBFixer cell to fix assets/1cgl_no_lig.pdb (adds missing hydrogens). Outputs assets/1cgl_fixed.pdb.  
Step 2: Run subsequent cells to perform MD simulations (define vacuum ChemicalSystem for ligands, set simulation settings [1 ns production run], create Protocol & NonTransformation, and write JSON files to run in terminal)  
- NOTE: system is performed in vacuum to improve runtime. In creating the settings, the settings for running NVT and NPT are purposely set to be really low to reduce runtime (our MD study only concerned about RMSD).
Step 3. Run RMSD cell to generate a plot of Backbone RMSD vs. Time, along with average and max RMSD values for each ligand  

## Expected Outputs

For Molecular Docking:
- `results/log_ligand{i}.txt` (AutoDock Vina output showing binding affinity scores, 9 poses per ligand)
- `results/output_ligand_{i}.pdbqt` (docked poses for ligands)
- `results/output_visual_ligand_{i}.sdf` (docked sdf poses)
- `table_results/table_ligand_{i}.csv` (table showing affinity, rmsd_lb, and rmsd_ub)
- Bar chart showing binding affinity of tested ligands (kcal/mol)
- RMSD upper bound line chart graphs RMSD across all 9 poses for the 3 candidate ligands (and the control)  

For Molecular Dynamics:
- MAIN OUTPUT `rmsd.png` Backbone RMSD vs. Time plot
- Mean and max RMSD values for all ligands (found in last cell of notebook)
- `md_results` (each ligand folder contains: equil_npt.pdb, equil_nvt.pdb, minimized.pdb, simulation.log)

## Runtime Notes And Limitations

For MD, each ligand takes ~30-60 minutes on CPU (for a 1ns simulation). For our purposes, a 1 ns runtime is mostly sufficient (RMSD results plateus and mostly stabilizes). Though a GPU is not required to run, it would greatly speed up the runtime of the simulations. A vacuum environment is created to reduce runtime, and so RMSD values may not accurately reflect how the ligand wound perform in physiological solvent conditions. 

