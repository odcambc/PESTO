import MDAnalysis as mda
import argparse

# Set up the argument parser
parser = argparse.ArgumentParser(
    description="This script calculates a pseudodensity of a selection of atoms in a PDB file."
)
parser.add_argument(
    "-i",
    "--input",
    type=str,
    help="The input PDB file to calculate a pseudodensity from",
)
parser.add_argument(
    "-o",
    "--output",
    type=str,
    help="The name of the output DX file to save the pseudodensity to",
)
parser.add_argument(
    "-d",
    "--delta",
    type=float,
    help="The delta value to use for the density calculation (default: 0.5)",
    default=0.5,
)
cli_args = parser.parse_args()

input_file = cli_args.input
output_file = cli_args.output
delta = cli_args.delta

structure_only = mda.Universe(input_file)
bb = structure_only.select_atoms("protein and backbone")

x = mda.analysis.density.DensityAnalysis(bb, delta=delta)
x.run()
x.density.convert_density(unit="Angstrom^{-3}")
x.density.export(output_file, type="dx")


def align_density(sel, distance):
    """
    @param sel: str, selection of the atoms to create an alignment density for
    @param distance: float, the distance in angstroms to consider for the alignment density
    """

    stored.sel = []
    # full lists
    stored.mol = []
    sel = sel + " and N. CA"

    # Iterate over residues in the selection
    cmd.iterate(sel, "stored.sel.append((model, resi))")

    # Iterate over models
    models = cmd.get_object_list()
    for model in models:
        if model != stored.mol:
            cmd.iterate(model, "stored.mol.append((model, resi))")

    # Count the number of residues within a given distance of the
    # residues in the selection and return an array containing the counts for each residue in the selection
    # by iterating over the stored selection, selecting residues using pymol's "within" selection operator,
    # and counting the number of residues in the selection
    count = []
    for model, resi in stored.sel:
        cmd.select("tmp", f"{model} and resi {resi}")
        cmd.select("tmp2", f"({sel}) within {distance} of tmp")
        count.append(cmd.count_atoms("tmp2"))
        cmd.delete("tmp2")
        cmd.delete("tmp")

    return count


def aligned_density(sel, step_distance=1):
    """
    From a selection, iterates over the x,y,z coordinates of a box around the selection
    and, for each step, counts the number of CA atoms in all other models within half the step distance
    @param sel: str, selection of the atoms to create an alignment density for
    @param step_distance: float, the distance (in angstroms) each step should be
    """

    # Get the bounding box of the selection
    cmd.select("tmp", sel)
    bbox = cmd.get_extent("tmp")
    cmd.delete("tmp")

    # Get the number of steps in each direction
    x_steps = int((bbox[1][0] - bbox[0][0]) / step_distance)
    y_steps = int((bbox[1][1] - bbox[0][1]) / step_distance)
    z_steps = int((bbox[1][2] - bbox[0][2]) / step_distance)

    # Create a dictionary to store the counts
    counts = {}

    # Iterate over the bounding box
    for x in range(x_steps):
        for y in range(y_steps):
            for z in range(z_steps):
                # Get the coordinates of the current step
                x_coord = bbox[0][0] + x * step_distance
                y_coord = bbox[0][1] + y * step_distance
                z_coord = bbox[0][2] + z * step_distance

                # Create a selection for the current step
                cmd.select(
                    "tmp",
                    f"({sel}) and (x > {x_coord} and x < {x_coord + step_distance}) and (y > {y_coord} and y < {y_coord + step_distance}) and (z > {z_coord} and z < {z_coord + step_distance})",
                )
                if cmd.count_atoms("tmp") == 0:
                    count = 0
                else:
                    # Count the number of CA atoms in all other models within half the step distance
                    count = 0
                    for model in cmd.get_object_list():
                        cmd.select("tmp2", f"{sel} within {step_distance / 2} of tmp")
                        count += cmd.count_atoms("tmp2")
                        cmd.delete("tmp2")

                # Store the count in the dictionary
                counts[(x_coord, y_coord, z_coord)] = count

                # Delete the temporary selection
                cmd.delete("tmp")

        return counts


def aligned_density(sel, step_distance=1, file=""):
    """
    From a selection, iterates over the x,y,z coordinates of a box around the selection
    and, for each step, counts the number of CA atoms in all other models within half the step distance
    @param sel: str, selection of the atoms to create an alignment density for
    @param step_distance: float, the distance (in angstroms) each step should be
    """

    # Get the bounding box of the selection
    cmd.select("tmp", sel)
    bbox = cmd.get_extent("tmp")
    cmd.delete("tmp")

    # Get the number of steps in each direction
    x_steps = int((bbox[1][0] - bbox[0][0]) / step_distance)
    y_steps = int((bbox[1][1] - bbox[0][1]) / step_distance)
    z_steps = int((bbox[1][2] - bbox[0][2]) / step_distance)

    radius = step_distance / 2

    # Create a dictionary to store the counts
    counts = {}

    # Iterate over the bounding box
    for x in range(x_steps):
        for y in range(y_steps):
            for z in range(z_steps):
                # Get the coordinates of the current step
                x_coord = bbox[0][0] + x * step_distance
                y_coord = bbox[0][1] + y * step_distance
                z_coord = bbox[0][2] + z * step_distance

                # Create a selection for the current step
                cmd.select(
                    "tmp",
                    f"({sel}) and (x > {x_coord - radius} and x < {x_coord + radius}) and (y > {y_coord - radius} and y < {y_coord + radius}) and (z > {z_coord - radius} and z < {z_coord + radius})",
                )
                try:
                    atom_count = cmd.count_atoms("tmp")
                except NameError:
                    atom_count = 0
                    print(f"{x_coord}, {y_coord}, {z_coord}")
                except:
                    atom_count = 0
                    print("Error counting atoms")

                if atom_count == 0:
                    count = 0
                else:
                    if sel == "all":
                        count = cmd.count_atoms("tmp")
                    else:
                        count = cmd.count_atoms("tmp") - 1

                # Store the count in the dictionary
                counts[(x_coord, y_coord, z_coord)] = count

                # Delete the temporary selection
                cmd.delete("tmp")

        # Print to file
        if file != "":
            with open(file, "w") as f:
                print(counts, f)

        return counts


def get_density(step_distance=1):
    """
    From a selection, iterates over the x,y,z coordinates of a box around all models
    and, for each step, counts the number of CA atoms in all other models within half the step distance
    @param step_distance: float, the distance (in angstroms) each step should be
    """

    return aligned_density(all, step_distance=step_distance)
