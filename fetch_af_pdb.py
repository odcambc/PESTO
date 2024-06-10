import argparse
import sys
import subprocess
import csv

# read in args from command line using argparse
parser = argparse.ArgumentParser(
    description="This script will read a csv file consisting of AlphaFold identifiers and download them to a specified folder."
)
parser.add_argument("-i", "--input", type=str, help="The input identifier csv file")
parser.add_argument(
    "-o",
    "--output",
    type=str,
    help="The name of the output directory to save the PDB files to",
)
parser.add_argument(
    "-v",
    "--verbose",
    action="store_true",
    help="Verbose output (prints curl output to stdout)",
)
parser.add_argument(
    "-s",
    "--silent",
    action="store_true",
    help="Silent output (suppresses all output to stdout - stderr will still be printed if an error occurs)",
)
args = parser.parse_args()

input_file = args.input
output_dir = args.output
if args.verbose:
    curl_args = ""
else:
    curl_args = "-s"

i = 0

if not args.silent:
    with open(input_file, "rb") as f:
        n_pdbs = sum(1 for _ in f)
    print(f"Total number of PDB files to download: {n_pdbs}")
    print(f"Downloading PDB files to {output_dir}")


with open(input_file, "r", encoding="utf-8") as f:
    reader = csv.reader(f)
    for row in reader:
        alphafold_ID = row[0]
        model_url = f"https://alphafold.ebi.ac.uk/files/{alphafold_ID}.pdb"
        try:
            retcode = subprocess.call(
                f"curl {curl_args} {model_url} -o {output_dir}/{alphafold_ID}.pdb",
                shell=True,
            )
            if retcode < 0:
                print("Curl was terminated by signal", -retcode, file=sys.stderr)
            else:
                if not args.silent:
                    print(f"Downloaded {alphafold_ID}.pdb: {i} of {n_pdbs}")
        except OSError as e:
            print("Execution failed:", e, file=sys.stderr)
        i += 1

if not args.silent:
    print("Download complete!")
