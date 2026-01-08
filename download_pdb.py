def download_pdb_file(pdb_id, output_dir="."):
    """Download PDB file"""
    print("\n=== Downloading PDB file ===")

    # Ensure directory exists 
    os.makedirs(protein_directory, exist_ok=True)
    pdb_file = f"{protein_directory}/{pdb_id}.pdb"
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    # Create parameter text file
    with open("params.json", "w") as param_file:
        param_file.write('{"pdb_id": "' + pdb_id + '"}\n')

    try:
        urllib.request.urlretrieve(url, pdb_file)
        print(f"Downloaded: {pdb_file}")
        return pdb_file
    except Exception as e:
        print(f"Download error: {e}")
        exit()
