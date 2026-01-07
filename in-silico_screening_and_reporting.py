#!/usr/bin/env python
# coding: utf-8

# # Step 1 Installation of dependencies
# Gnina will run within a linux environment provided by google colab virtual machine.
# 
# 1. `useful_rdkit_utils` is a Python package written and maintained by Pat Walters that contains useful RDKit functions. We will use it for the functions `mcs_rmsd` (explained later).
# 2. `py3Dmol` is used for molecular visualization.
# 3. The RDKit is a popular cheminiformatics package we will use for processing molecules.


import subprocess
import os
import requests
import useful_rdkit_utils as uru
from rdkit import Chem
from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.rdBase import BlockLogs
import pandas as pd

### Check for previous gnina installation ###
gnina_file = "gnina"

# Download URL
gnina_url = "https://github.com/gnina/gnina/releases/download/v1.3/gnina"

# Check if file exists
if os.path.exists(gnina_file):
    print(f"✅ {gnina_file} already exists, skipping download.")
else:
    print(f"⬇️ {gnina_file} not found, downloading...")
    subprocess.run([
        "wget",
        gnina_url,
        "-O", gnina_file
    ], check=True)
    
    # Make gnina executable
    subprocess.run(["chmod", "+x", "gnina"])
    print("✅ Download complete.")

pdb_id = input("Enter PDB code used in protein_preparation.py: ")
ligand_id = input("Enter ligand code used in ligand_extraction_and_preparation.py: ")

### Working directories ###
protein_directory = "molecular_docking/protein_files"
ligand_directory = "molecular_docking/ligand_structures"
docking_results_directory = "molecular_docking/docking_results"
docking_results = ""
# ------------------------------------------------------------------------------
# 4.1 Docking modes
# ------------------------------------------------------------------------------

def single_ligand_docking():
    subprocess.run([
    "./gnina",
    "-r", f"{protein_directory}/{pdb_id}_A.pdbqt",
    "-l", f"{ligand_directory}/{ligand_id}.sdf",
    "--autobox_ligand", f"{ligand_directory}/{ligand_id}_corrected_pose.sdf",
    "-o", f"{docking_results_directory}/{ligand_id}_docked_{pdb_id}.sdf",
    "--seed", "0",
    "--exhaustiveness", f"{ex}",
    *gpu_flag,  # expands to nothing or ["--no_gpu"]
    *cnn_flag  # expands to nothing or ["--cnn_scoring=none"]
    ])

def batch_docking():
    subprocess.run([
    "./gnina",
    "-r", f"{protein_directory}/{pdb_id}_A.pdbqt",
    "-l", f"{ligand_directory}/ligands_to_dock.sdf",
    "--autobox_ligand", f"{ligand_directory}/{ligand_id}_corrected_pose.sdf",
    "-o", f"{docking_results_directory}/multiple_ligands_docked_{pdb_id}.sdf",
    "--seed", "0",
    "--exhaustiveness", f"{ex}",
    *gpu_flag,  # expands to nothing or ["--no_gpu"]
    *cnn_flag  # expands to nothing or ["--cnn_scoring=none"]
    ])

def flexible_docking():
    subprocess.run([
    "./gnina",
    "-r", f"{protein_directory}/{pdb_id}_A.pdbqt",
    "-l", f"{ligand_directory}/{ligand_id}.sdf",
    "--autobox_ligand", f"{ligand_directory}/{ligand_id}_corrected_pose.sdf",
    "-o", f"{docking_results_directory}/{ligand_id}_flex.sdf",
    "--flexdist_ligand", f"{ligand_directory}/{ligand_id}_corrected_pose.sdf",
    "--flexdist", "3.59",
    "--seed", "0",
    "--exhaustiveness 64"
    *gpu_flag,  # expands to nothing or ["--no_gpu"]
    *cnn_flag  # expands to nothing or ["--cnn_scoring=none"]
    ])

def unknown_site_docking():
    subprocess.run([
    "./gnina",
    "-r", f"{protein_directory}/{pdb_id}_A.pdbqt",
    "-l", f"{ligand_directory}/{ligand_id}.sdf",
    "--autobox_ligand", f"{protein_directory}/{pdb_id}_A.pdbqt",
    "-o", f"{docking_results_directory}/{ligand_id}_docked_whole_{pdb_id}.sdf",
    "--seed", "0",
    "--exhaustiveness", f"{ex}",
    *gpu_flag,  # expands to nothing or ["--no_gpu"]
    *cnn_flag  # expands to nothing or ["--cnn_scoring=none"]
    ])
    
# ------------------------------------------------------------------------------
# 4.2 Main functions
# ------------------------------------------------------------------------------

def docking_main():
    # Redocking with extracted ligand
    if selection == "a":
        print("\n === Single ligand docking ===")
        single_ligand_docking()
        rmsd_calculation(ligand_directory, docking_results_directory, ligand_id, pdb_id)
        report(docking_results = f"{docking_results_directory}/{ligand_id}_docked_{pdb_id}.sdf")
    
    # Docking with multiple ligands
    elif selection == "b":
        print("\n === Batch docking ===")
        batch_docking()
        report(docking_results = f"{docking_results_directory}/multiple_ligands_docked_{pdb_id}.sdf")
    
    # Flexible docking
    elif selection == "c":
        print("Flexible docking")
        flexible_docking()
        report(docking_results = f"{docking_results_directory}/{ligand_id}_flex.sdf")
    
    # Docking on unknown site
    else:
        print("Docking on unknown site")
        unknown_site_docking()
        report(docking_results = f"{docking_results_directory}/{ligand_id}_docked_whole_{pdb_id}.sdf")
        
def report(docking_results):
    score_columns = [
        "minimizedAffinity",
        "CNNscore",
        "CNNaffinity",
        "CNN_VS",
        "CNNaffinity_variance",
    ]

    # Normalize input: accept either a string or a list of strings
    if isinstance(docking_results, str):
        sdf_paths = [docking_results]
    else:
        sdf_paths = docking_results
    
    df_list = []
    for filename in sdf_paths:
        with BlockLogs():
            df_list.append(PandasTools.LoadSDF(filename))
    
    combo_df = pd.concat(df_list)
    
    # Convert score columns to float if they exist and are not entirely empty
    for col in score_columns:
        if col in combo_df.columns:
            # If the column is completely empty, drop it
            if combo_df[col].isnull().all() or (combo_df[col] == "").all():
                combo_df.drop(columns=[col], inplace=True)
            else:
                combo_df[col] = combo_df[col].astype(float)
    
    # Extract SMILES from the RDKit molecule column
    combo_df["SMILES"] = combo_df["ROMol"].apply(lambda mol: Chem.MolToSmiles(mol) if mol is not None else None)
    
    # Save with SMILES included
    output_csv = f"{docking_results_directory}/docking_results_{pdb_id}.csv"
    combo_df.to_csv(output_csv, index=False)
    
    print(f"Docking results with SMILES saved to {output_csv}")
    
    combo_df.head()

# ------------------------------------------------------------------------------
# 4.3 Additional functions
# ------------------------------------------------------------------------------

def rmsd_calculation(ligand_directory, docking_results_directory, ligand_id, pdb_id):
    print("\n === Running mcs (maximum common structure) rmsd calculation ===")
    cognate = Chem.MolFromMolFile(f"{ligand_directory}/{ligand_id}_corrected_pose.sdf")
    poses = Chem.SDMolSupplier(f"{docking_results_directory}/{ligand_id}_docked_{pdb_id}.sdf")

    # Prepare log file path
    log_path = os.path.join(docking_results_directory, f"{ligand_id}_{pdb_id}_rmsd.log")

    with open(log_path, "w") as log_file:
        log_file.write("Pose_Index\tNum_Matches\tRMSD\n")  # header line

        for i, pose in enumerate(poses):
            if pose is None:
                continue
            RDLogger.DisableLog('rdApp.warning')
            n_match, rmsd = uru.mcs_rmsd(cognate, pose)
            line = f"{i}\t{n_match}\t{rmsd:.2f}\n"
            print(line.strip())       # print to console
            log_file.write(line)      # save to log file

    print(f"\n === RMSD results saved to {log_path} ===")

# Run gnina
import os

os.makedirs("molecular_docking/docking_results", exist_ok=True)

selection = input("\n === Welcome to gnina. Please select the docking mode: \
                    \n single docking (rmsd and cnn score will be calculated): a \
                    \n batch docking: b \
                    \n flexible docking (maximum exhaustiveness): c \
                    \n docking on unknown sites: d \
                ")

inp = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80]

ex = int(input("\n Define exhaustiveness (8, 16, 24, 32, 40, 48, 56, 64): "))
if ex in inp:
    print(f"You have selected  exhaustiveness level {ex}")
    ask_gpu = input("\n Run with GPU? [y/n]: ").strip().lower()
    ask_cnn = input("\n Run with CNN? [y/n]: ").strip().lower()
    
    # Decide GPU flag
    gpu_flag = []
    if ask_gpu == 'y':
        gpu_flag = []  # run with GPU
    elif ask_gpu == 'n':
        gpu_flag = ["--no_gpu"]  # run without GPU
    else:
        print("\n Invalid input. Please try again with 'y' or 'n'.")
    
    # Decide CNN flag
    cnn_flag = []
    if ask_cnn == 'y':
        cnn_flag = []  # run with CNN
    elif ask_cnn == 'n':
        cnn_flag = ["--cnn_scoring=none"]  # run without CNN
    else:
        print("\n Invalid input. Please try again with 'y' or 'n'.")

    print(f"You have selected option {selection}.")

    if __name__ == "__main__":
        #pdb_id = input("Enter PDB code used in protein_preparation.py: ")
        #pdb_id = os.getenv("PARAM_PDB_ID")
        #ligand_id = input("Enter ligand code used in ligand_extraction_and_preparation.py: ")
        #ligand_name = os.getenv("PARAM_LIGAND_ID")

        #with open("params.json") as f:
            #data = json.load(f)

        #ligand_id = data["ligand_id"]
        #print("Ligand ID parsed:", ligand_id)

        docking_main()
else:
    print("\n === Invalid input. Please try again with the given options only! ===")
    print("\n === Exiting gnina ===")
    ex = ''




