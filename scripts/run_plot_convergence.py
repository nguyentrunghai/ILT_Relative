"""
plot convergence of pearson's R and RMSE with respect to number of receptor snapshots
"""

from __future__ import print_function

import argparse
import os

import numpy as np

from load_mbar_weights_holo_OBC2 import load_mbar_weights
from _plots import improved_plot_lines

parser = argparse.ArgumentParser()

parser.add_argument("--data_1_dir", type=str, default="Relative_FE_Est_1")
parser.add_argument("--data_2_dir", type=str, default="Relative_FE_Est_with_CV_method_3a")
parser.add_argument("--data_file", type=str, default="r_rmse.dat")

# "pearson_R" or "RMSE"
parser.add_argument("--which_data", type=str, default="none")

# to modify data_2
parser.add_argument("--ref_ligands_2_modify", type=str, default=" ")
parser.add_argument("--modifying_constants", type=str, default=" ")

parser.add_argument("--xlabel", type=str, default="# receptor snapshots")
parser.add_argument("--ylabel", type=str, default="Pearson's R w.r.t. final results")

parser.add_argument("--colors", type=str, default="g b")
parser.add_argument("--line_styles", type=str, default="- --")
parser.add_argument("--line_width", type=float, default=2)

args = parser.parse_args()


def _load_data(data_file, which_data):
    data = np.loadtxt(data_file)
    if which_data == "pearson_R":
        size, r, r_std = data[:, 0],  data[:, 1], data[:, 2]
        return size, r, r_std
    elif which_data == "RMSE":
        size, rmse, rmse_std = data[:, 0],  data[:, 3], data[:, 4]
        return size, rmse, rmse_std
    else:
        raise ValueError("Unknown which_data: " + which_data)


ref_ligands_2_modify = args.ref_ligands_2_modify.split()
modifying_constants = [np.float(s) for s in args.modifying_constants.split()]
assert len(ref_ligands_2_modify) == len(modifying_constants), "len(ref_ligands_2_modify) neq len(modifying_constants)"

colors = args.colors.split()
line_styles = args.line_styles.split()
assert len(colors) == len(line_styles) == 2, "len(colors) and len(line_styles) must equal 2"

modifying_constants = {ligand : num for ligand, num in zip(ref_ligands_2_modify, modifying_constants)}
print("modifying_constants", modifying_constants)

_, _, single_snap_weights, _, _ = load_mbar_weights()

ref_ligands = [ligand for ligand in single_snap_weights.keys() if ligand != "systems"]
for ref_ligand in ref_ligands:
    data_file_1 = os.path.join(args.data_1_dir, ref_ligand, args.data_file)
    print("Loading " + args.which_data + " from " + data_file_1)
    data_1 = _load_data(data_file_1, args.which_data)

    data_file_2 = os.path.join(args.data_2_dir, ref_ligand, args.data_file)
    print("Loading " + args.which_data + " from " + data_file_2)
    data_2 = _load_data(data_file_2, args.which_data)

    xs = [data_1[0], data_2[0]]

    y1 = data_1[1]
    y2 = data_2[1]
    if ref_ligand in modifying_constants:
        y2 += modifying_constants[ref_ligand]
    ys = [y1, y2]

    yerr1 = data_1[2] / 2.
    yerr2 = data_2[2] / 2.
    yerrs = [yerr1, yerr2]

    out = ref_ligand + "_" + args.which_data + ".pdf"
    print("Plotting "+out)
    improved_plot_lines(xs, ys, yerrs=yerrs, xlabel=args.xlabel, ylabel=args.ylabel, out=out,
                        colors=colors,
                        line_styles=line_styles,
                        lw=args.line_width)
