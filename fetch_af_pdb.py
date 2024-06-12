import argparse
import sys
from pathlib import Path
import requests


def test_url(model_url):
    """Test if the formatted URL is valid. Curl will get a 404 error code
    if the pdb file does not exist, otherwise 200.
    """
    test_response = requests.get(model_url, timeout=5)
    if test_response.ok:
        return True
    else:
        print(f"Error: request status code {test_response.status_code}")
        return False


def fetch_model(model_url):
    """Fetch a model from the AlphaFold database using the provided URL."""
    model_request = requests.get(model_url, timeout=5)
    if model_request.ok:
        return model_request.content
    else:
        return None


def check_file(model_id, dir):
    """Check if the file already exists in the output directory."""
    if Path(f"{dir}/{model_id}.pdb").is_file():
        return True
    else:
        return False


def write_request(model_request, dir, model_id):
    """Write the model request content to a file."""
    with open(f"{dir}/{model_id}.pdb", "wb") as output_file:
        output_file.write(model_request)


# read in args from command line using argparse
parser = argparse.ArgumentParser(
    description="This script will read a FoldSeek output hit file in m8 format, extract the identifiers, and download them to a specified folder."
)
parser.add_argument("-i", "--input", type=str, help="The input m8 file")
parser.add_argument("--format", type=str, help="Format of input file (m8 or list)")
parser.add_argument(
    "-o",
    "--output",
    type=str,
    help="The name of the output directory to save the PDB files to",
)
parser.add_argument(
    "-s",
    "--silent",
    action="store_true",
    help="Silent output (no progress messages)",
)
parser.add_argument(
    "-f",
    "--force_download",
    action="store_true",
    help="Force download of all files (overwrite existing files)",
)
parser.add_argument(
    "--min_qlen",
    type=int,
    help="Skip models with qlen below this value",
)
parser.add_argument(
    "--min_tlen",
    type=int,
    help="Skip models with tlen below this value",
)
args = parser.parse_args()

input_file = args.input

# Set input format
if args.format:
    INPUT_FORMAT = args.format
else:
    if input_file.endswith(".m8"):
        INPUT_FORMAT = "m8"
    elif input_file.endswith(".csv"):
        INPUT_FORMAT = "list"
    elif input_file.endswith(".txt"):
        INPUT_FORMAT = "list"
    else:
        print("Input file format not recognized. Please specify format with --format")
        sys.exit(1)
if not args.output:
    output_dir = Path.cwd()
else:
    output_dir = Path(args.output)

i = 0

if not args.silent:
    with open(input_file, "rb") as f:
        n_pdbs = sum(1 for _ in f)
    print(f"Total number of PDB files to download: {n_pdbs}")
    print(f"Downloading PDB files to {output_dir}")


with open(input_file, "r", encoding="utf-8") as f:
    # test if parsing works before downloading any more files
    start_pos = f.tell()
    test_line = f.readline()
    if INPUT_FORMAT == "list":
        test_alphafold_ID = test_line.strip()
    elif INPUT_FORMAT == "m8":
        test_alphafold_ID = test_line.split()[1]
    url = f"https://alphafold.ebi.ac.uk/files/{test_alphafold_ID}.pdb"
    test_url_response = test_url(url)
    if not test_url_response:
        sys.exit(1)
    else:
        if not args.silent:
            print(f"Tested URL {url} successfully.")
    f.seek(start_pos)

    # Fetch the rest of the models
    for line in f:
        if INPUT_FORMAT == "list":
            alphafold_ID = line.strip()
        elif INPUT_FORMAT == "m8":
            alphafold_ID = line.split()[1]
            # Check for qlen and tlen if specified
            if args.min_qlen:
                if int(line.split('\t')[9]) < args.min_qlen:
                    if not args.silent:
                        print(f"Skipping {alphafold_ID} due to low sequence identity (qlen).")
                    continue
            if args.min_tlen:
                if int(line.split('\t')[12]) < args.min_tlen:
                    if not args.silent:
                        print(f"Skipping {alphafold_ID} due to low sequence identity (tlen).")
                    continue
        if not args.force_download:
            if check_file(alphafold_ID, output_dir):
                if not args.silent:
                    print(f"File {alphafold_ID}.pdb already exists. Skipping.")
                i += 1
                continue
        else:
            if not args.silent:
                print(f"Redownloading {alphafold_ID}.pdb")

        url = f"https://alphafold.ebi.ac.uk/files/{alphafold_ID}.pdb"
        request = fetch_model(url)
        if request:
            write_request(request, output_dir, alphafold_ID)
            if not args.silent:
                print(f"Downloaded {alphafold_ID}.pdb: {i} of {n_pdbs}")
        else:
            # If the request is not ok, print the error message but don't stop.
            print(
                f"Error: request status code {request.status_code} for ID {alphafold_ID}."
            )
        i += 1

if not args.silent:
    print("Download complete!")
