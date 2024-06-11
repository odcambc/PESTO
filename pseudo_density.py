import argparse
import MDAnalysis as mda
import MDAnalysis.analysis.density

# read in args from command line using argparse
parser = argparse.ArgumentParser(
    description="This script will take a a pdb file and generate a pseduodensity map in dx format."
)
parser.add_argument("-i", "--input", type=str, help="The input pdb file")
parser.add_argument("-o", "--output", type=str, help="The name of the output dx file")
parser.add_argument(
    "-d",
    "--delta",
    type=float,
    help="The grid spacing for sampling PDB density (in Angstroms), default is 0.5",
    default=0.5,
)
parser.add_argument(
    "--ca", "--CA", action="store_true", help="Use CA (backbone) atoms only"
)
args = parser.parse_args()

input_file = args.input
output_file = args.output
delta = args.delta

pdb_input = mda.Universe(input_file)
if args.ca:
    pdb_atoms = pdb_input.select_atoms("protein and backbone")
else:
    pdb_atoms = pdb_input.select_atoms("protein")

x = MDAnalysis.analysis.density.DensityAnalysis(pdb_atoms, delta=delta)
x.run()
x.density.convert_density(unit="Angstrom^{-3}")
x.density.export(output_file)
