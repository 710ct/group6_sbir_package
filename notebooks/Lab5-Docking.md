# Lab 5: Molecular Docking with AutoDock Vina
### University of California, Berkeley - Spring 2026 - ME 120/292A
---
## Overview

### Purpose
The purpose of this lab is to run a small-molecule docking workflow end-to-end using AutoDock Vina, from downloading structures, receptor + ligand preparation, docking, and visualization of poses and binding scores.

### Environment
In this lab, we will use a modern, purely command-line (CLI) workflow for maximum reliability across different operating systems:
- **AutoDock Vina** (https://vina.scripps.edu/) for docking and scoring ligand poses.
- **PDBFixer** (https://github.com/openmm/pdbfixer) to clean raw protein structures and add missing hydrogen atoms.
- **Meeko** (https://github.com/forlilab/Meeko) for converting both ligand and receptor structures into the required `.pdbqt` format.
- **Mol*** (https://molstar.org/viewer/) for modern, web-based 3D visualization of the receptor, grid box, and final docking poses to capture publication-quality screenshots.
- **AutoDockTools** (https://ccsb.scripps.edu/mgltools/) as an optional, legacy graphical interface for visualizing docking parameters and preparing structures (demonstrated in supplementary material).


*(Note: The traditional graphical tool, AutoDockTools, is notoriously buggy on modern systems and DataHub. It will be demonstrated in a supplementary video for context, but you are not required to install or use it for this lab).*

Files you will work with:
- Terminal, command-line, text editor (`.txt`, `.pdb`, `.sdf`, `.pdbqt`)

> Platforms: UC Berkeley DataHub or local (macOS, Linux, Windows (WSL)).

### Objectives
By the end of the lab, you should be able to:
- Download a receptor structure from the Protein Data Bank and a ligand structure from PubChem
- Convert and prepare ligand and receptor structures into `.pdbqt` format
- Define a docking search box around a binding site
- Run AutoDock Vina using a configuration file and interpret docking scores
- Visualize docking poses and capture publication-quality screenshots

### Submission
Submit a single PDF report to bCourses `lab5_{firstname}_{lastname}.pdf`.
* **Part A (4 pts)**: Guided Example Analysis (beta2AR + albuterol).
* **Part B (6 pts)**: Independent Bacterial System Analysis.
(See the Deliverable Assignment section at the end of the document for exact requirements).

---
## Setup

### Files
Download the lab directory from the bCourses **Labs/** folder as a `.zip` archive. 

**Extract** this archive into your `ME120/` or `ME292A/` course folder on your DataHub or local filesystem. You should see the following directory:
- `Lab5-Docking/`

Navigate into this folder in your terminal before starting the lab.

### Dependencies
For this lab, we will use two separate Conda environments. The **Required** environment (`vina`) contains all the modern command-line tools you need to complete the assignment. The **Optional** environment (`mgltools`) contains the legacy AutoDockTools graphical interface, which your instructor will use in the supplementary video. 

*(Note: MGLTools relies on older software architecture, which is why it must be kept in its own isolated environment to prevent conflicts).*

#### 1. The Main Workflow Environment (Required)
- **Conda env name:** `vina`

**Create the environment from YAML (Recommended):**
```bash
conda env create -f lab5-environment.yml
```

Create environment manually (Alternative):
```bash
conda create -n vina python=3.10 -y
conda activate vina
# Add conda-forge channel if you have not already
conda config --env --add channels conda-forge
conda install numpy scipy rdkit swig boost-cpp libboost sphinx sphinx_rtd_theme gemmi autogrid vina meeko openmm pdbfixer -y
pip install vina prody
```

Activate the environment:
```bash
conda activate vina
```
> Keep this environment activated for the entirety of the docking assignment.

#### 2. (Optional) AutoDockTools GUI Environment 
- Conda env name: `mgltools`

If you wish to follow along with the video and use the traditional graphical interface to visualize your grid boxes, create this secondary environment:
  
Create the environment manually
```bash
conda create -n mgltools -y
conda activate mgltools
conda install bioconda::mgltools -y
```

Start AutoDockTools from the terminal:
```bash
adt
```

> Note: AutoDockTools requires a working display/GUI. It may not launch properly on remote headless servers like DataHub without specific X11 forwarding configurations

#### Switching between environments

```bash
conda activate vina      # required 
conda activate mgltools  # optional GUI
```

---
## 1. Introduction

Molecular docking is a computational technique used to predict the preferred orientation and binding affinity of a small molecule (ligand) when it interacts with a protein (receptor). This process helps simulate the molecular interactions between the ligand and the receptor, such as hydrogen bonding, hydrophobic interactions, and electrostatic forces. The goal of docking is to identify how well a ligand fits into the receptor’s binding site and predict the most energetically favorable binding mode. Docking is widely used in drug discovery to screen potential therapeutic compounds, optimize binding interactions, and design molecules with higher specificity and efficacy for their target proteins.

<p align="center">
  <img src="https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/examples.jpg" alt="AutoDock Vina Examples" width="500">
</p>

## 2. Docking beta2AR with albuterol

### Background

In this assignment, we will focus on the Beta-2 Adrenergic Receptor (β2AR), a protein found on the surface of cells in the lungs, heart, and other tissues. This receptor plays a key role in the regulation of smooth muscle relaxation, particularly in the airways. **Salbutamol**, referred to as **Albuterol** in North America and sold under proprietary names such as Ventolin, Proventil, or ProAir, is a drug commonly used to treat asthma, a chronic respiratory disease that causes inflammation and constriction of the airways, leading to difficulty breathing.

By binding to β2AR, Salbutamol/Albuterol causes bronchodilation, or the widening of the airways, which helps alleviate asthma symptoms and improve airflow. Using molecular docking, we can explore how Salbutamol/Albuterol interacts with the β2AR at the molecular level, providing insights into the drug’s binding mechanism and its therapeutic effects.

This process is essential in the development of new and more effective asthma medications, helping researchers design drugs that can better target the receptor and optimize therapeutic outcomes. In this lab, we will learn how to acquire receptor and ligand structures, prepare them for docking using **PDBFixer**, for protein cleanup, **Meeko** for **PDBQT** preparation, run docking with **AutoDock Vina**, and analyze the resulting poses.

Optionally, in the accompanying video, we also demonstrate **AutoDockTools** as a legacy GUI option for visualizing the binding site and understanding docking box parameters, but it is not required for the assignment.

<table style="width:100%; text-align:center; border:none;">
  <tr>
    <td><img src="https://cdn.rcsb.org/images/structures/3p0g_model-1.jpeg" alt="Protein" style="max-height:250px;"></td>
    <td><img src="https://upload.wikimedia.org/wikipedia/commons/thumb/d/d5/RS-salbutamol-from-xtal-3D-balls.png/250px-RS-salbutamol-from-xtal-3D-balls.png" alt="Ligand" style="max-height:250px;"></td>
    <td><img src="https://s.turbifycdn.com/aah/yhst-135855760451349/ventolin-hfa-albuterol-sulfate-90mcg-200-metered-inhalations-54.jpg" alt="Inhaler" style="max-height:250px;"></td>
  </tr>
</table>

In Part 1 you will dock albuterol (salbutamol) to the beta-2 adrenergic receptor (beta2AR).

### Download Receptor and Ligand Files

#### Receptor: beta2AR (PDB ID: 3P0G)
1. Go to the [RCSB Protein Data Bank](https://www.rcsb.org/).
2. Search for PDB ID: `3P0G`.
3. Download the structure in **Legacy PDB** format.
4. Save the file into your lab folder as `beta2ar.pdb`.

*Optional: view the known ligand binding site on RCSB*
1. On the `3P0G` page, find the **Biological Assembly 1** panel.
2. Click **Ligand Interaction** to visualize where the crystallographic ligand binds.
3. Use this visualization to guide your docking search box placement in the deliverable.

#### Ligand: albuterol (PubChem CID: 2083)
1. Go to [PubChem](https://pubchem.ncbi.nlm.nih.gov/).
2. Search for PubChem CID: `2083`.
3. Download the ligand in **SDF** format.
4. Save the file into your lab folder as `albuterol.sdf`.

> Checkpoint: you have `beta2ar.pdb` and `albuterol.sdf` in working directory.

*Optional: download both files from the command line directly*
```bash
curl -L -o beta2ar.pdb https://files.rcsb.org/download/3P0G.pdb
curl -L -o albuterol.sdf "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/2083/SDF?record_type=3d"
```

<div style="display:flex; flex-wrap:wrap; gap:12px; align-items:stretch;">

  <!-- Left: Protein -->
  <div style="flex:1 1 420px;">
    <div style="border:1px solid #ddd; border-radius:8px; overflow:hidden;">
      <iframe
        src="https://molstar.org/viewer/?pdb=3P0G&pdb-provider=rcsb&hide-controls=1&collapse-left-panel=1"
        width="100%"
        height="520"
        style="border:0;"
        loading="lazy"
        title="Mol* Preview: 3P0G (beta2AR)">
      </iframe>
    </div>
    <p style="margin:8px 0 0 0;">
      <a href="https://molstar.org/viewer/?pdb=3P0G&pdb-provider=rcsb" target="_blank" rel="noopener">
        Open beta2AR (3P0G) in Mol* (new tab)
      </a>
    </p>
  </div>

  <!-- Right: Ligand from PubChem (3D SDF) -->
  <div style="flex:1 1 420px;">
    <div style="border:1px solid #ddd; border-radius:8px; overflow:hidden;">
      <iframe
        src="https://molstar.org/viewer/?hide-controls=1&collapse-left-panel=1&structure-url=https%3A%2F%2Fpubchem.ncbi.nlm.nih.gov%2Frest%2Fpug%2Fcompound%2FCID%2F2083%2FSDF%3Frecord_type%3D3d&structure-url-format=sdf"
        width="100%"
        height="520"
        style="border:0;"
        loading="lazy"
        title="Mol* Preview: Albuterol (PubChem CID 2083)">
      </iframe>
    </div>
    <p style="margin:8px 0 0 0;">
      <a href="https://molstar.org/viewer/?structure-url=https%3A%2F%2Fpubchem.ncbi.nlm.nih.gov%2Frest%2Fpug%2Fcompound%2FCID%2F2083%2FSDF%3Frecord_type%3D3d&structure-url-format=sdf"
         target="_blank" rel="noopener">
        Open albuterol (CID 2083) in Mol* (new tab)
      </a>
    </p>
  </div>

</div>

### Prepare the ligand (`ligand.pdbqt`)

AutoDock Vina requires input files to be in the `.pdbqt` format, which includes partial charges and atom types. We will use **Meeko**, the official Python package developed by the AutoDock suite creators, to automatically add hydrogens, calculate charges, and convert our SDF file.
1. Open a **terminal** and ensure your `vina` conda environment is activated.
2. Run the Meeko preparation script on your SDF file:
```bash
mk_prepare_ligand.py -i albuterol.sdf -o ligand.pdbqt
```
> Checkpoint: The file `ligand.pdbqt` has been generated in your working directory.


### Prepare the receptor (`receptor.pdbqt`)

Proteins downloaded from the RCSB PDB often lack hydrogen atoms and contain unwanted elements like water molecules or crystallization artifacts.

First, we will use PDBFixer to clean the structure and add missing hydrogens at a physiological pH (7.4):
```bash
pdbfixer beta2ar.pdb \
  --output beta2ar_receptorH.pdb \
  --keep-heterogens=none \
  --replace-nonstandard \
  --ph=7.4
```

Next, we will use Meeko to assign the correct AutoDock atom types, create the .pdbqt file, and simultaneously generate our docking search box configuration file:


```bash
mk_prepare_receptor.py -i beta2ar_receptorH.pdb -o receptor -p -v \
  --box_center 64.610 17.949 12.734 \
  --box_size   25.0   25.0   25.0
```

- `-i beta2ar_receptorH.pdb`: input receptor structure (cleaned).

- `-o receptor`: output basename (creates receptor.*).

- `-p`: write `receptor.pdbqt` (the final docking receptor).

- `-v`: write `receptor.box.txt` (the Vina box config).

- `--box_center`, `--box_size`: search box center and dimensions (in Ångströms).

> Checkpoint: `receptor.pdbqt` and `receptor.box.txt` have been created in your working directory.


### Run docking with AutoDock Vina 

Before running docking, it is useful to see what options Vina provides.
```bash
vina --help | less
```
Key options:
- `--receptor`: rigid part of the protein receptor (PDBQT).

- `--ligand`: small-molecule ligand (PDBQT).

- `--center_x/y/z` and `--size_x/y/z`: dimensions of the search box, in Å.

- `--exhaustiveness`: exhaustiveness of the global search, roughly proportional to time.

- `--seed`: random seed for repeatability.

- `--out`: output models (PDBQT), the default is chosen based on the ligand file name

To run a docking simulation, AutoDock Vina needs several pieces of information: what the protein is, what the drug is, where to look, and how hard to search.

If we wanted to run Vina entirely through the command line, we would have to use a long list of flags:
```bash
vina --receptor receptor.pdbqt --ligand ligand.pdbqt \
  --config receptor.box.txt \
  --exhaustiveness 9 --seed 100 \
  --out output.pdbqt
```

### Best Practice: Using a `config.txt` file

For reproducibility, it is better to store all run settings in a single `config.txt` file. Then you can re-run the same “experiment” later with one command, or make a controlled change by copying the config and editing one line.

Create a text file named `config.txt` in your working directory with the following content.
Replace `{center_x}`, `{center_y}`, `{center_z}`, `{size_x}`, `{size_y}`, `{size_z}` with your box parameters generated by **Meeko**.

```text
receptor = receptor.pdbqt
ligand = ligand.pdbqt

# Center of the search box (Angstrom)
center_x = {center_x}
center_y = {center_y}
center_z = {center_z}

# Box size (Angstrom)
size_x = {size_x}
size_y = {size_y}
size_z = {size_z}

# Increase for longer / more "accurate" docking (slower)
exhaustiveness = 9

# Set a seed for repeatability
seed = 100
```

Run docking using only the config file:
```bash
vina --config config.txt --out output.pdbqt
```

Since all of your parameters are stored in one file, your command becomes simple.

One final note: Vina prints the final affinity scores directly to your terminal screen (stdout). Newer versions of Vina do not have a built-in `--log` option to save this text into a file automatically. Therefore, we must use standard terminal commands to capture the output ourselves into a log file.

You have two choices for how to run your simulation:

1. Use the overwrite redirect `>` to push to a log file:
```bash
vina --config config.txt --out output.pdbqt > log.txt
```

2. (Recommended) Use the pipe `|` operator to feed the stdout into the `tee` command. Just like a "T" junction, `tee` splits the output: it prints it to your terminal screen so you can watch the progress bar, and it saves it to a file at the same time.
```bash
vina --config config.txt --out output.pdbqt | tee log.txt
```
> Checkpoint: `output.pdbqt` and `log.txt` were created.

### Analyze and View Results

Before visualizing the 3D structure, look at the numeric results.

1. Open `log.txt` and find the table of docking modes.
2. Record the binding affinities (kcal/mol) for the best modes. Mode 1 is the most energetically favorable (lowest / most negative affinity).

To view how the ligand actually fits into the protein pocket, you have two options:

#### Option 1: Mol* Viewer (Robust)
Mol* (Molstar) can load PDB/PDBQT and SDF, and the [Mol* viewer](molstar.org/viewer/) has an option that lets you load multiple structures into the same scene.
1. Open your web browser and go to: **[https://molstar.org/viewer/](https://molstar.org/viewer/)**
2. On the left side, click **Open Files** and upload the following three files from your lab folder simultaneously:
   * `beta2ar_receptorH.pdb` (cleaned protein structure)
   * `receptor.box.pdb` (The 3D grid box generated by Meeko)
   * `output.pdbqt` (The 9 docking poses from Vina)
3. **Visualize the Poses:** By default, Mol* loads all 9 ligand poses on top of each other. 
   * On the right-hand control panel, locate the `output.pdbqt` component.
   * Look for the **Trajectory** or **Model** slider controls. Use the arrows to cycle through models 1 through 9 to see how the different poses of the ligand.
4. **Style your Screenshot:** For a publication-quality render:
   * Set the `beta2ar_receptorH.pdb` representation to **Cartoon** (to show the protein ribbons).
   * Set the `output.pdbqt` representation to **Ball & Stick** or **Sticks**.
5. Position your camera to clearly show the ligand inside the binding pocket and the bounding box, then take a screenshot of the best mode (1/9) for your report using the circular aperture button in the right side of the view.

#### Option 2 (Optional): AutoDockTools
If you created the optional `mgltools` conda environment and are running the lab locally, you can use the classic interface:

1. Launch AutoDockTools from your terminal by typing `adt` in your `mgltools` env.
2. **Load the protein:** Go to **File** > **Read Molecule** and open your `receptor.pdbqt`.
3. Click the **R** (oval button) option on the dashboard next to the protein name to switch to ribbon view.
4. **Load the poses:** Open your `output.pdbqt` file.
5. **Cycle the poses:** Enable the conformations tool in ADT and use the arrow keys to cycle through the 9 different poses to identify the best mode.
6. Take a screenshot of the best mode as your deliverable.



## 3. Flexible Docking

Rigid docking treats the receptor as perfectly static, but real binding sites are not rigid. Side chains near the pocket can rotate, shift slightly to relieve steric clashes, or form/break hydrogen bonds as a ligand settles into place. **AutoDock Vina** supports limited receptor flexibility by allowing one or more side chains to move during docking. Flexible docking is a small step towards a more realistic model; it keeps the overall protein rigid for speed, but allows key side chains to move to see whether a plausible local rearrangement changes the predicted pose or affinity.

In this part, we will repeat our docking simulation but treat one specific binding-site residue as flexible.
### Prepare a flexible receptor

This command splits the receptor into two PDBQT files:

* `flex_receptor_rigid.pdbqt` (most of the receptor)
* `flex_receptor_flex.pdbqt` (the flexible side chain)

We will treat **Ser203 (chain A, residue 203)** as flexible.

**Why Serine 203?** Biologically, Ser203 sits right in the middle of the binding pocket. When an asthma drug like albuterol binds, it forms a strong hydrogen bond with this specific Serine, physically pulling it slightly to "activate" the receptor. Because we know this residue moves in real life, making it flexible in our simulation allows Vina to find a better hydrogen-bonding angle. 

In short, this information was obtained from literature. For your independent bacterial system in Part B, you will need to look at primary literature or the RCSB PDB entry to identify a logical side chain to make flexible. 

```bash
mk_prepare_receptor.py -i beta2ar_receptorH.pdb -o flex_receptor -p -v \
  --box_center 64.610 17.949 12.734 \
  --box_size   25.0   25.0   25.0 \
  -f A:203 -a
```

* `-f A:203` selects the flexible residue (chain A, residue 203).
* `-a` is the optional ignore/remove residues that fail template matching safeguard

> Checkpoint: you have both `flex_receptor_rigid.pdbqt` and `flex_receptor_flex.pdbqt`, plus `flex_receptor.box.txt`.

### Run flexible docking
> This uses the vina forcefield

Run Vina using your existing `config.txt` file, but use the command line flags to override the receptor with your new **rigid** part and the **flex** side chain:

```bash
vina --config config.txt \
  --receptor flex_receptor_rigid.pdbqt \
  --flex flex_receptor_flex.pdbqt \
  --exhaustiveness 9 \
  --out output_flex.pdbqt | tee log_flex.txt
```

> Checkpoint: `output_flex.pdbqt` and `log_flex.txt` were created.

### Compare to rigid docking (quick check)

Open `log.txt` (rigid) and `log_flex.txt` (flexible) and compare the **best affinity (mode 1)**. You will include both results in your report.


### (FIX) Convert Flexible Results to SDF
Newer versions of the Mol* viewer hosted on their website may have inconsistencies with how the bond connectivity is determined. This issue seems most common with the **Flexible** results. We can use Meeko to convert our poses to the common `.sdf` format, which can hold multiple poses in one file. 

```bash
mk_export.py output_flex.pdbqt -s output_flex_poses_withflex.sdf -k
```

---
## Deliverable Assignment (10 pts)

Submit a short, typeset PDF report (**1–2 pages**) to bCourses named:

* `lab5_{firstname}_{lastname}.pdf`

Your report must include the required figures with captions and brief written responses. Refer to figures by number (Figure 1, Figure 2, etc.) in your text.

---
### Bacterial System Requirement (Part B)

For the independent portion of this lab, you must analyze a **receptor from a bacterial species** (domain: Bacteria) and dock a **small molecule ligand** to it.


#### System Choice Suggestions:
- **Source organism**: Bacteria (e.g., E. coli, S. aureus, M. tuberculosis)
- **Receptor**: choose a protein structure available for download from the RCSB PDB. Keep it under **800 residues** for ease.
- **Ligand**: choose a small molecule available for download from PubChem (or similar) in SDF format. Keep it under **500 Da (g/mol)** for ease. 
> It is highly recommended to choose a PDB structure that already has a co-crystallized drug or ligand bound to it. This makes defining your search box much easier!

### Required Figures
Include these figures in your PDF (use Mol* or AutoDockTools for your screenshots):
- **Figure 1**: Screenshot of the best rigid pose for the guided example (beta2AR + albuterol).
- **Figure 2**: Screenshot of the best flexible pose for the guided example.
- **Figure 3**: Screenshot of the best rigid pose for your chosen bacterial system.
- **Figure 4**: Screenshot of the best flexible pose for your chosen bacterial system.

### Part A (4 pts): Guided Example Analysis
*Analyze the beta2AR + albuterol docking results from the tutorial.*

#### Q1 (4 pts): Comparison and Interpretation - Rigid vs. Flexible

Write a brief response (**75–300 words**) addressing the following:
- Did allowing **Ser203** to be flexible change the best affinity? (Which is more favorable / more negative?)
- Did the pose/orientation of albuterol change noticeably between rigid and flexible docking?
- Why might local flexibility be important for simulating this specific drug-receptor interaction? *Hint: Think about "induced fit"*
- Include one short summary of whether the rest of the modes look similar or more variable between rigid vs flexible (you do not need to analyze all 9 in detail).

Reference **Figure 1** and **Figure 2** in your report. Provide the best affinity score (mode 1) for both the rigid run and the flexible run.

### Part B (6 pts): Independent Bacterial System

Choose **one bacterial receptor-ligand system** and repeat the docking workflow from this lab (download --> prepare ligand/receptor --> define simulation parameters --> run Vina --> analyze results).

**System requirements**

* Receptor: a **bacterial** protein structure (PDB ID required).
* Ligand: a small molecule ligand (PubChem CID or equivalent source).
* Your docking box must be justified (how you chose `center_*` and `size_*`).

#### Q2 (2 pts): System choice + box justification

Provide the setup details for the system you researched. Justify your choices for doing so. (**50-300 words**):

* Receptor name + **PDB ID** + organism (bacterial species)
* Ligand name + **PubChem CID** (or other source)
* A justification of how you chose the docking box and the parameters used:
  * `center_x`, `center_y`, `center_z`
  * `size_x`, `size_y`, `size_z`
* Which specific residue did you choose to make flexible and why?

#### Q3 (4 pts): Comparison and Interpretation - Rigid vs. Flexible - (Your Bacterial System)

Write a brief response (**200-400 words**) addressing the following:
- What is the normal function of this bacterial protein, and why and how is the ligand relevant? (Is it a known antibiotic? An inhibitor? A metabolite? etc.)
- Report the best affinity (Mode 1) for your rigid and flexible docking. Did adding flexibility change the binding score or the physical pose for your specific system? Try to provide an explanation on why it did or did not have an effect.
- Did the pose change meaningfully between rigid and flexible? Describe what changed by describing the location in pocket, orientation, contacts with residues.

Include at least one literature citation (a primary research paper, the PDB citation, or a reputable database page).

Reference **Figure 3** and **Figure 4** in your report. Provide the best affinity score (mode 1) for both the rigid run and the flexible run. You should also comment on any differences with the rest of the 9 modes for both simulations. 

---
## Troubleshooting

### Quick checks
If you run into issues, try these first:
1. Restart the terminal, then run again.
2. Confirm your conda environment is active: `conda env list`.
3. Confirm you are in the correct directory: `pwd` and `ls`.
4. Re-run the last command and read the first error message carefully.

### Common issues
- Issue: `adt` launches but the window is blank or crashes.
  Fix: try the separate `mgltools` conda env, or install MGLTools natively for your OS.

### Common issues (files and formats)
- Issue: Vina says it cannot read `receptor.pdbqt` or `ligand.pdbqt`.
  Fix: confirm the files exist and the names match `config.txt`. Re-run Meeko or re-export from AutoDockTools.

- Issue: Docking results look nonsensical (ligand docks far away).
  Fix: your search box is likely misplaced or too large. Re-check the known binding pocket and tighten the box.

- Issue: `vina: command not found`.
  Fix: ensure the conda env is activated (`conda activate vina`) and re-run the install command.
---
## Collaboration and Academic Integrity

You may discuss high-level ideas and debugging strategies with classmates.
Do not share completed reports, finished `.pdbqt` outputs, or written answers.
All submitted work must be your own.

Cite any external resources used (links are fine).
If you received help, include a short acknowledgment in your report.

---
## References

Lab links used in this document:
- AutoDock Vina - docking engine and manual - https://vina.scripps.edu/ and https://vina.scripps.edu/manual/
- MGLTools / AutoDockTools - structure preparation and grid box selection - https://ccsb.scripps.edu/mgltools/
- RCSB Protein Data Bank - receptor structures - https://www.rcsb.org/
- PubChem - ligand structures - https://pubchem.ncbi.nlm.nih.gov/
