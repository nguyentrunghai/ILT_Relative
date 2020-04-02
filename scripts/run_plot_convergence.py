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

parser.add_argument("--ref_systems",
type=str, default="1-methylpyrrole.A__AAA benzene.A__AAA lysozyme.active.A__ABJ lysozyme.inactive.A__AAS phenol.A__AAA p-xylene.A__AAA")

# "Pearson_R" or "RMSE"
parser.add_argument("--which_data", type=str, default="none")

parser.add_argument("--y2_scale_facs", type=str, default="")
parser.add_argument("--y2_shifts", type=str, default="")
parser.add_argument("--y2err_scale_facs", type=str, default="")

parser.add_argument("--xlabel", type=str, default="# receptor snapshots")
parser.add_argument("--ylabel", type=str, default="Pearson's R w.r.t. final results")

parser.add_argument("--colors", type=str, default="g b")
parser.add_argument("--line_styles", type=str, default="- --")
parser.add_argument("--line_width", type=float, default=2)

args = parser.parse_args()


def _load_data(data_file, which_data):
    data = np.loadtxt(data_file)
    if which_data == "Pearson_R":
        size, r, r_std = data[:, 0],  data[:, 1], data[:, 2]
        return size, r, r_std
    elif which_data == "RMSE":
        size, rmse, rmse_std = data[:, 0],  data[:, 3], data[:, 4]
        return size, rmse, rmse_std
    else:
        raise ValueError("Unknown which_data: " + which_data)


def scaling(x, fac=1.):
    x = x[::-1]
    facs = np.linspace(1, fac, len(x))
    x = x*facs
    return x[::-1]


def plot_convergence(data_file_1, data_file_2, which_data, ax,
                     y2_scale=1, y2err_scale=1, y2_shift=0,
                     colors=("b", "r"), line_styles=["-", "--"],
                     lw=1):
    x1, y1, yerr1 = _load_data(data_file_1, which_data)
    x2, y2, yerr2 = _load_data(data_file_2, which_data)

    yerr1 = yerr1 / 2
    yerr2 = yerr2 / 2
    yerr2 = yerr2 * y2err_scale

    y2 = scaling(y2, fac=y2_scale) + y2_shift

    ax.errorbar(x1, y1, yerr=yerr1, c=colors[0], linestyle=line_styles[0], lw=lw)
    ax.errorbar(x2, y2, yerr=yerr2, c=colors[1], linestyle=line_styles[1], lw=lw)

    return None


colors = args.colors.split()
line_styles = args.line_styles.split()
assert len(colors) == len(line_styles) == 2, "len(colors) and len(line_styles) must equal 2"

ref_systems = args.ref_systems.split()
print("ref_systems:", ref_systems)

# y2_scale
y2_scale_facs = args.y2_scale_facs.split()
if len(y2_scale_facs) == 0:
    y2_scale_facs = [1. for _ in ref_systems]
else:
    assert len(y2_scale_facs) == len(ref_systems), "wrong len of y2_scale_facs"
    y2_scale_facs = [float(s) for s in y2_scale_facs]
print("y2_scale_facs")

# y2_shifts
y2_shifts = args.split()
if len(y2_shifts) == 0:
    y2_shifts = [0. for _ in ref_systems]
else:
    assert len(y2_shifts) == len(ref_systems), "wrong len of y2_shifts"
    y2_shifts = [float(s) for s in y2_shifts]
print("y2_shifts:", y2_shifts)

# y2err_scale_facs
y2err_scale_facs = args.y2err_scale_facs.split()
if len(y2err_scale_facs) == 0:
    y2err_scale_facs = [1. for _ in ref_systems]
else:
    assert len(y2err_scale_facs) == len(ref_systems), "wrong len of y2err_scale_facs"
    y2err_scale_facs = [float(s) for s in y2err_scale_facs]
print("y2err_scale_facs:", y2err_scale_facs)



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
