#!/bin/bash

echo "=== Starting workflow ==="

# Download PDB
echo "Running download_pdb.py..."
python download_pdb.py
if [ $? -ne 0 ]; then
    echo "Error running download_pdb.py. Exiting!"
    exit 1
fi

# Protein preparation
echo "Running protein_preparation.py..."
python protein_preparation.py
if [ $? -ne 0 ]; then
    echo "Error running protein_preparation.py. Exiting!"
    exit 1
fi

# Ligand preparation
echo "Running ligand_extraction_and_preparation.py..."
python ligand_extraction_and_preparation.py
if [ $? -ne 0 ]; then
    echo "Error running ligand_extraction_and_preparation.py. Exiting!"
    exit 1
fi

# Docking
echo "Running in-silico_screening_and_reporting.py..."
python in-silico_screening_and_reporting.py
if [ $? -ne 0 ]; then
    echo "Error running in-silico_screening_and_reporting.py. Exiting!"
    exit 1
fi

echo "=== All scripts finished ==="
