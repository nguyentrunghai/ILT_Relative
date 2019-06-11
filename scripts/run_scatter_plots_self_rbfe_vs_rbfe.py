"""
run scatter plots between self relative binding free energies and relative binding free energies
TODO set figure titles
"""
from __future__ import print_function

import os
import glob
import argparse

import numpy as np

from _plots import scatter_plot

parser = argparse.ArgumentParser()

parser.add_argument("--data_dir", type=str,
                    default="/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/test/p-xylene.A__AAA__equal_sys__single_weight/ExpMean")

parser.add_argument("--glob_matching", type=str, default="*_vs_*")

parser.add_argument("--xlabel", type=str, default="Self relative binding free energy (kcal/mol)")
parser.add_argument("--ylabel", type=str, default="Relative binding free energy (kcal/mol)")

args = parser.parse_args()

data_files = glob.glob(os.path.join(args.data_dir, args.glob_matching))

for data_file in data_files:
    print("Processing file", data_file)
    data = np.loadtxt(data_file)
    x = data[:, 0]
    y = data[:, 1]

    out_file = os.path.basename(data_file) + ".pdf"

    scatter_plot(x, y, args.xlabel, args.ylabel, out_file,
                 show_xy_axes=True,
                 xerr=None, yerr=None,
                 show_regression_line=True,
                 show_diagonal_line=False,
                 show_rmse=True,
                 show_R=True,
                 show_regression_line_eq=True,
                 markers=None,
                 markersize=4,
                 markercolors=None,
                 same_xy_scale=False,
                 text_pos=[0.1, 0.7],
                 title=os.path.basename(data_file))
