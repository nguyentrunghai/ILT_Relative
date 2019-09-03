"""
plot convergence of pearson's R and RMSE with respect to number of receptor snapshots
"""

from __future__ import print_function

import argparse
import os

import nump as np

from _yank import YANK_LIGANDS as ref_ligands

parser = argparse.ArgumentParser()

parser.add_argument("--data_1_dir", type=str, default="Relative_FE_Est_1")
parser.add_argument("--data_2_dir", type=str, default="Relative_FE_Est_with_CV_method_3a")
parser.add_argument("--data_file", type=str, default="r_rmse.dat")

# "pearson_R" or "RMSE"
parser.add_argument("--which_data", type=str, default="none")


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


for ref_ligand in ref_ligands:
    data_file_1 = os.path.join(args.data_1_dir, ref_ligand, args.data_file)
    print("Loading " + args.which_data + " from " + data_file_1)
    data_1 = _load_data(data_file_1, args.which_data)

    data_file_2 = os.path.join(args.data_2_dir, ref_ligand, args.data_file)
    print("Loading " + args.which_data + " from " + data_file_2)
    data_2 = _load_data(data_file_2, args.which_data)
    