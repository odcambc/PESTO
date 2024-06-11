import argparse
import subprocess
import itertools
import threading
import queue
from pathlib import Path
import pandas as pd

# read in args from command line using argparse
parser = argparse.ArgumentParser(
    description="This script perform pairwise alignment of all pdb files in a folder using usalign and save TM-scores."
)
parser.add_argument(
    "-i", "--indir", type=str, help="The input directory containing PDB files"
)
parser.add_argument(
    "-o",
    "--output",
    type=str,
    help="The name of the output csv file to save the results to",
)
parser.add_argument(
    "-d",
    "--aligned_dir",
    type=str,
    help="The location of the directory to save aligned pdb files to (default: ./aligned)",
)
parser.add_argument(
    "-j",
    "--threads",
    type=int,
    help="The number of threads to use (default: 8)",
    default=8,
)
parser.add_argument(
    "-r",
    "--reference",
    type=str,
    help="The reference structure to align all other structures to",
)
parser.add_argument(
    "-q",
    "--quiet",
    action="store_true",
    help="Suppress output progress (default: False)",
)
parser.add_argument(
    "-s",
    "--silent",
    action="store_true",
    help="Suppress all output, but still show errors (default: False)",
)
parser.add_argument(
    "-x",
    "--executable",
    type=str,
    help="Specify the path to the USalign executable (default: ~/bin/USalign)",
)
cli_args = parser.parse_args()

if cli_args.silent:
    QUIET = True
    SILENT = True
elif cli_args.quiet:
    QUIET = True
    SILENT = False
else:
    QUIET = False
    SILENT = False

# Check if directory is supplied. If not, use current directory.
if not cli_args.indir:
    DIR = "."
else:
    DIR = cli_args.indir

# Check if directory exists.
if not Path(DIR).is_dir():
    if not SILENT:
        print(f"Error: directory {DIR} not found")
        exit(1)
else:
    if not SILENT:
        print(f"Aligning files in directory {DIR}")

if cli_args.executable:
    USALIGN = Path(cli_args.executable).expanduser().resolve()
    # Check if USalign runs
    try:
        subprocess.call(USALIGN, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        print(f"{USALIGN} not found")
        exit(1)
else:
    try:
        subprocess.call("USalign", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        USALIGN = "USalign"
    except FileNotFoundError:
        # handle file not found error.
        print("USalign not found in path")
        try:
            USALIGN = Path("~/bin/USalign").expanduser().resolve()
            subprocess.call(str(USALIGN), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except FileNotFoundError:
            print("USalign not found at ~/bin/USalign")
            print("Please specify the path to the USalign executable using the -x flag")
            exit(1)

OUT_FILE = cli_args.output
NUM_WORKER_THREADS = cli_args.threads
aligned_output_dir = cli_args.aligned_dir

reference_model = cli_args.reference
# Check reference file: if it doesn't end with pdb, add it and then check.
if reference_model:
    # Check reference model name
    if not reference_model.endswith(".pdb"):
        reference_model = reference_model + ".pdb"
    if not Path(f"{DIR}/{reference_model}").is_file():
        print(f"{DIR}/{reference_model}")
        if not SILENT:
            print(f"Reference model {reference_model} not found in {DIR}")
        exit(1)
    if not SILENT:
        print(f"Using reference model {reference_model}")
    # Check aligned output folder
    if not aligned_output_dir:
        aligned_output_dir = Path("./aligned")
    else:
        aligned_output_dir = Path(aligned_output_dir)
    if not aligned_output_dir.is_dir():
        print(aligned_output_dir)
        aligned_output_dir.mkdir(parents=True)
        if not SILENT:
            print(f"Created directory {str(aligned_output_dir)}")
    else:
        aligned_output_dir = Path(aligned_output_dir)
        if not SILENT:
            print(f"Saving aligned files to {str(aligned_output_dir)}")


def parse_usalign_output(stdout):
    structure1 = ""
    length1 = 0
    structure2 = ""
    length2 = 0
    aligned_length = 0
    rmsd = 0
    seq_id = 0
    tmscore_n1 = 0
    tmscore_n2 = 0
    alignment1 = ""
    alignment2 = ""
    pairs = ""

    line = str(next(stdout))

    while not str(line).startswith("b'Name of Structure_1:"):
        line = next(stdout)

    # Need to parse both of these type of outputs
    # 'Name of Structure_1: AF-Q58918-F1-model_v4.pdb:A (to be superimposed onto Structure_2)\n'
    # 'Name of Structure_1: 5FXB/AF-Q32IV5-F1-model_v4.pdb:A (to be superimposed onto Structure_2)\n'
    # This is jank but works
    structure1 = str(line).split(":")[1].split("/")[-1].split(".")[0].strip()
    structure2 = str(next(stdout)).split(":")[1].split("/")[-1].split(".")[0].strip()
    length1 = int(str(next(stdout)).split(":")[1].split(" ")[1])
    length2 = int(str(next(stdout)).split(":")[1].split(" ")[1])
    next(stdout)
    line = str(next(stdout)).split(" ")
    aligned_length = int(line[2][:-1])
    rmsd = float(line[6][:-1])
    seq_id = float(line[8][:-3])
    tmscore_n1 = float(str(next(stdout)).split(" ")[1])
    tmscore_n2 = float(str(next(stdout)).split(" ")[1])
    next(stdout)
    next(stdout)
    next(stdout)
    alignment1 = str(next(stdout))
    pairs = str(next(stdout))
    alignment2 = str(next(stdout))

    parsed_output = {
        "structure1": structure1,
        "length1": length1,
        "structure2": structure2,
        "length2": length2,
        "aligned_length": aligned_length,
        "rmsd": rmsd,
        "seq_id": seq_id,
        "tmscore_n1": tmscore_n1,
        "tmscore_n2": tmscore_n2,
        "alignment1": alignment1,
        "pairs": pairs,
        "alignment2": alignment2,
    }

    return parsed_output


def run_usalign_and_parse(pdb_pair):
    # Runs USalign and just generates the TM-score
    command = f"{str(USALIGN)} {pdb_pair[0]} {pdb_pair[1]}"
    cmd = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    parsed_output = parse_usalign_output(cmd.stdout)

    try:
        tm_scores_df.loc[parsed_output["structure1"], parsed_output["structure2"]] = (
            parsed_output["tmscore_n1"]
        )
    except KeyError:
        print(
            f"Error with {parsed_output['structure1']} and {parsed_output['structure2']}"
        )
        print(parsed_output)


def run_usalign_and_parse_with_aligned(pdb_pair):
    # Runs USalign with option "-o" to output an aligned pdb file
    # -o superimposes structure1 onto structure2
    aligned_pdb_name = pdb_pair[0].stem
    reference_pdb_name = pdb_pair[1].stem
    mapped_pdb_name = f"{aligned_pdb_name}_onto_{reference_pdb_name}"

    output_file = aligned_output_dir / mapped_pdb_name

    command = f"{str(USALIGN)} {pdb_pair[0]} {pdb_pair[1]} -o {str(output_file)}"
    cmd = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    parsed_output = parse_usalign_output(cmd.stdout)

    try:
        tm_scores_df.loc[parsed_output["structure1"], parsed_output["structure2"]] = (
            parsed_output["tmscore_n1"]
        )
    except KeyError:
        print(
            f"Error with {parsed_output['structure1']} and {parsed_output['structure2']}"
        )
        print(parsed_output)


pdb_file_list = [x for x in Path(DIR).iterdir() if x.suffix == ".pdb"]
n_files = len(pdb_file_list)

pdb_name_list = [pdb.stem for pdb in pdb_file_list]
combinations = list(itertools.combinations(pdb_file_list, 2))
n_comparisons = len(combinations)
tm_scores_df = pd.DataFrame(columns=pdb_name_list, index=pdb_name_list)


class Worker(threading.Thread):
    def __init__(self, qq, *args, **kwargs):
        self.q = qq
        super().__init__(*args, **kwargs)

    def run(self):
        while True:
            try:
                pair = self.q.get(timeout=3)  # 3s timeout
            except queue.Empty:
                return
            if pair[1].name == f"{reference_model}":
                run_usalign_and_parse_with_aligned(pair)
            else:
                run_usalign_and_parse(pair)
            if not QUIET:
                print(f"{self.q.qsize()} / {n_comparisons} remaining...")
            self.q.task_done()


q = queue.Queue()
for combination in combinations:
    q.put_nowait(combination)
for _ in range(NUM_WORKER_THREADS):
    Worker(q).start()
q.join()  # blocks until the queue is empty.

tm_scores_df.to_csv(OUT_FILE)
if not SILENT:
    print("done!")
