import argparse
import subprocess
import itertools
import threading
import queue
import os
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
    "-j",
    "--threads",
    type=int,
    help="The number of threads to use (default: 8)",
    default=8,
)
cli_args = parser.parse_args()


DIR = cli_args.indir
OUT_FILE = cli_args.output
NUM_WORKER_THREADS = cli_args.threads


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

    structure1 = str(line).split(":")[1].split("/")[1].split(".")[0]
    structure2 = str(next(stdout)).split(":")[1].split("/")[1].split(".")[0]
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
    command = f"~/bin/USalign {pdb_pair[0]} {pdb_pair[1]}"
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


pdb_file_list = [f"{DIR}/" + x for x in os.listdir(DIR) if x.endswith(".pdb")]

pdb_name_list = [pdb.split("/")[1].split(".")[0] for pdb in pdb_file_list]
combinations = list(itertools.combinations(pdb_file_list, 2))
total = len(combinations)
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
            run_usalign_and_parse(pair)
            print(f"{self.q.qsize()} / {total}")
            self.q.task_done()


q = queue.Queue()
for combination in combinations:
    q.put_nowait(combination)
for _ in range(NUM_WORKER_THREADS):
    Worker(q).start()
q.join()  # blocks until the queue is empty.

tm_scores_df.to_csv(OUT_FILE)
print("done!")
