import argparse
from pathlib import Path
import subprocess

# read in args from command line using argparse
parser = argparse.ArgumentParser(
    description="This script combines all PDB files in a directory into a single one."
)
parser.add_argument(
    "-i", "--indir", type=str, help="The input directory containing PDB files to merge"
)
parser.add_argument(
    "-o",
    "--output",
    type=str,
    help="The name of the output PDB file to save the results to",
)
parser.add_argument(
    "-r",
    "--reference",
    type=str,
    help="The reference structure all structures are aligned to",
)
cli_args = parser.parse_args()

if not cli_args.indir:
    print("Please specify an input directory with the -i flag.")
    exit(1)
if not cli_args.output:
    print("Please specify an output file with the -o flag.")
    exit(1)
if not cli_args.reference:
    print("Please specify the reference PDB file with the -r flag.")
    exit(1)

input_dir = "aligned"
input_dir = "."
aligned_file = "aligned_pdbs.pdb"
reference_pdb = "../5FXB/AF-O25823-F1-model_v4.pdb"

# Set up variables
input_dir = cli_args.indir
aligned_file = cli_args.output
reference_pdb = cli_args.reference


# Get list of PDBs to combine
pdb_file_list = [str(x) for x in Path(input_dir).iterdir() if x.suffix == ".pdb"]

# Tidy reference PDB and save as aligned file with pdb-tools
cmd = subprocess.run(
    f"pdb_tidy {reference_pdb} > {aligned_file}", shell=True, check=True
)

# Combine all PDBs into a single PDB
for pdb in pdb_file_list:
    # Copy the aligned file to a temp file, then merge
    cmd = subprocess.run(["cp", f"{aligned_file}", f"{aligned_file}.tmp"], shell=False)
    cmd = subprocess.run(
        f"pdb_merge {aligned_file}.tmp {pdb} > {aligned_file}", shell=True, check=True
    )

# Clean up temp file
cmd = subprocess.call(["rm", f"{aligned_file}.tmp"], shell=False)
