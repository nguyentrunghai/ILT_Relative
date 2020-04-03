"""
plot convergence of pearson's R and RMSE with respect to number of receptor snapshots
"""

from __future__ import print_function

import argparse
import os

import numpy as np

import matplotlib.pyplot as plt
#import seaborn as sns
#sns.set()

parser = argparse.ArgumentParser()

parser.add_argument("--data_1_dir", type=str, default="Relative_FE_Est_1")
parser.add_argument("--data_2_dir", type=str, default="Relative_FE_Est_with_CV_method_3a")
parser.add_argument("--data_file", type=str, default="r_rmse.dat")

parser.add_argument("--ref_systems",
type=str, default="1-methylpyrrole.A__AAA benzene.A__AAA lysozyme.active.A__ABJ lysozyme.inactive.A__AAS phenol.A__AAA p-xylene.A__AAA")

# "Pearson_R" or "RMSE"
parser.add_argument("--which_data", type=str, default="none")

parser.add_argument("--y2_scale_facs", type=str, default="")
parser.add_argument("--y1_shifts", type=str, default="")
parser.add_argument("--y2_shifts", type=str, default="")
parser.add_argument("--y2err_scale_facs", type=str, default="")

parser.add_argument("--xlabel", type=str, default="# receptor snapshots")
parser.add_argument("--ylabel", type=str, default="Pearson's R w.r.t. final results")

parser.add_argument("--colors", type=str, default="b r")
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
                     y2_scale=1., y2err_scale=1., y1_shift=0., y2_shift=0.,
                     colors=("b", "r"), line_styles=("-", "--"),
                     lw=1):
    x1, y1, yerr1 = _load_data(data_file_1, which_data)
    x2, y2, yerr2 = _load_data(data_file_2, which_data)

    yerr1 = yerr1 / 2
    yerr2 = yerr2 / 2
    yerr2 = yerr2 * y2err_scale

    y1 = y1 + y1_shift
    y2 = scaling(y2, fac=y2_scale) + y2_shift

    ax.errorbar(x1, y1, yerr=yerr1, c=colors[0], linestyle=line_styles[0], lw=lw)
    ax.errorbar(x2, y2, yerr=yerr2, c=colors[1], linestyle=line_styles[1], lw=lw)

    return None


colors = args.colors.split()
line_styles = args.line_styles.split()
assert len(colors) == len(line_styles) == 2, "len(colors) and len(line_styles) must equal 2"

line_width = args.line_width

ref_systems = args.ref_systems.split()
print("ref_systems:", ref_systems)

# y2_scale
y2_scale_facs = args.y2_scale_facs.split()
if len(y2_scale_facs) == 0:
    y2_scale_facs = [1. for _ in ref_systems]
else:
    assert len(y2_scale_facs) == len(ref_systems), "wrong len of y2_scale_facs"
    y2_scale_facs = [float(s) for s in y2_scale_facs]
print("y2_scale_facs:", y2_scale_facs)

# y1_shifts
y1_shifts = args.y1_shifts.split()
if len(y1_shifts) == 0:
    y1_shifts = [0. for _ in ref_systems]
else:
    assert len(y1_shifts) == len(ref_systems), "wrong len of y1_shifts"
    y1_shifts = [float(s) for s in y1_shifts]
print("y1_shifts:", y1_shifts)

# y2_shifts
y2_shifts = args.y2_shifts.split()
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


which_data = args.which_data
print("which_data:", which_data)

FONTSIZE = 8
FONT = {"fontname": "Arial"}

for i, ref_system in enumerate(ref_systems):
    data_file_1 = os.path.join(args.data_1_dir, ref_system, args.data_file)
    print("data_file_1:", data_file_1)

    data_file_2 = os.path.join(args.data_2_dir, ref_system, args.data_file)
    print("data_file_2:", data_file_2)

    y2_scale = y2_scale_facs[i]
    y1_shift = y1_shifts[i]
    y2_shift = y2_shifts[i]
    y2err_scale = y2err_scale_facs[i]

    fig, ax = plt.subplots(1, 1, figsize=(3.2, 2.4))
    plot_convergence(data_file_1, data_file_2, which_data, ax,
                     y2_scale=y2_scale, y2err_scale=y2err_scale,
                     y1_shift=y1_shift, y2_shift=y2_shift,
                     colors=colors, line_styles=line_styles, lw=line_width)

    ax.set_xlim([8, 100])

    ax.set_xlabel(args.xlabel, fontsize=FONTSIZE, **FONT)
    ax.set_ylabel(args.ylabel, fontsize=FONTSIZE, **FONT)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(FONTSIZE)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(FONTSIZE)

    fig.tight_layout()
    out = ref_system + "_" + which_data + ".pdf"
    print("Writing figure to " + out)
    fig.savefig(out)

print("DONE")
