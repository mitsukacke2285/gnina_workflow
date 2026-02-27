# Simplified Docking with Gnina (SDoG) - a fully automated Gnina workflow
This simple workflow for docking small molecule ligands into receptors is ideal for hit-to-lead or lead optimization. Its core modules are:

1. Download of protein crystal structure of interest
2. Protein preparation with PDBFixer, rdkit, MDAnalysis, Biopython
3. Ligand preparation with Molscrub
4. Docking with Gnina (based on Autodock Vina and Smina)

This workflow can be used by beginner and advanced users alike. The only input needed is the PDB id. The ligand(s) will be automatically identified after entering the PDB id in the ligand preparation module.

Scope and limitations:
1. rigid and flexible docking (later: free energy perturbation)
2. small molecule ligand docking
3. docking into known binding site of co-crystallized ligand
4. docking into unknown binding site(s)

# Requirements/packages needed to be installed
Biopython, Gnina, MDAnalysis, Numpy, OpenBabel-Wheel, OpenMM, Os, Pandas, PDBFixer, Pytest, Requests, Rdkit, Molscrub, Scipy, Subprocess, Useful-rdkit-utils, a list of SMILES of compounds to be prepared for docking as a csv file, Docker Desktop

PLEASE NOTE: This workflow was only tested in WSL (Windows Subsystem for linux)!!! There is no guarantee that it runs in other OS!!!

# Installation
Create the necessary environment for SDoG to run. In WSL type the following commands:

1. docker build -t {your_username}/workflow-gnina:full . ---> builds image from full.dockerfile
2. docker push {your_username}/workflow-gnina:full ---> pushes image to your Docker Hub account
3. docker run -it -v $(pwd):/workspace {your_username}/workflow-gnina:full bash ---> creates and runs container

# Usage
You can either run "bash runall.sh" to execute the entire workflow automatically (download_pdb.py, protein_preparation.py, ligand_extraction_and_preparation.py, in-silico_screening_and_reporting.py) or run each .py file separately. The following section describes the workflow of each py files:

## Download PDB file

download_pdb.py

Script will prompt you to enter your desired PDB code, downloads it and creates the file "params.json" which stores the PDB name used for the entire workflow. 

## Protein preparation

python protein_preparation.py

The working folder "molecular_docking/protein_files". The target protein will be prepared automatically and a pdbqt file will be generated at the end of the script.

## Ligand preparation

python ligand_extraction_and_preparation.py

The working folder "molecular_docking/ligand_structures" will be created. The script will identify all ligands bound to the PDB structure and prompts you to select a ligand you are interested in. This ligand will be saved as a ligand_id variable and will stick with the rest of the script. Several ligand files (sdf, pdb) will be generated along the way. The final output is:

1. a ligand that is pose corrected against an ideal ligand (downloaded automatically from https://www.rcsb.org/ and saved as {ligand_id}_ideal.sdf) and prepared by molscrub and saved as {ligand_id}.sdf.

2. an sdf file containing several scrubbed and pose corrected compounds obtained from the csv file.

3. "params.json" will be updated with the name of the extracted ligand for the remainder of the workflow.

## Docking with Gnina

python in-silico_screening_and_reporting.py

Gnina will be downloaded and installed, and a folder named molecular_docking/docking_results will be created. You will be prompted to select the docking mode by typing in the desired letter (selections a-d). After some time when the script finishes, the output will be a docked ligand sdf and a report csv file (name output depends on the docking mode selected earlier). The results can be analyzed in an external program or script.

# References
- PDB101 tutorial by RCSB Protein Data Bank: https://pdb101.rcsb.org/train/training-events/python4
- McNutt, A.T., Li, Y., Meli, R. et al. GNINA 1.3: the next increment in molecular docking with deep learning. J Cheminform 17, 28 (2025). https://doi.org/10.1186/s13321-025-00973-x
- Buccheri, R.; Rescifina, A. High-Throughput, High-Quality: Benchmarking GNINA and AutoDock Vina for Precision Virtual Screening Workflow. Molecules 2025, 30, 3361. https://doi.org/10.3390/molecules30163361
