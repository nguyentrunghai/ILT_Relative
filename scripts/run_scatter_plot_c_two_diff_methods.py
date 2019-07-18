"""
plot histograms of c for individual methods and scatter plots for pairs of methods
"""

from __future__ import print_function

import argparse
import glob
import os

import numpy as np

from _algdock import SIX_YANK_SYSTEMS

parser = argparse.ArgumentParser()
parser.add_argument("--data_dirs", type=str,
                    default="../Relative_FE_Est_with_CV_method_2a/not_flip_sign_c ../Relative_FE_Est_with_CV_method_2b/not_flip_sign_c")
parser.add_argument("--method_labels", type=str,
                    default="method_2a method_2b")
args = parser.parse_args()

data_dirs = args.data_dirs.split()
print("data_dirs", data_dirs)
method_labels = args.method_labels.split()
print("method_labels")

ref_ligands = SIX_YANK_SYSTEMS
SUB_DIR_SUFFIX = "__equal_sys__single_weight"

cs = {}
for label, data_dir in zip(method_labels, data_dirs):
    cs[label] = {}

    for ref_ligand in ref_ligands:
        cs[label][ref_ligand] = {}

        data_files = glob.glob(os.path.join(data_dir, ref_ligand + SUB_DIR_SUFFIX, "ExpMean",
                                            ref_ligand + "_G_CORR_H_*"))
        data_files = [data_file for data_file in data_files if data_file.split("_G_CORR_H_")[-1] != ref_ligand]

        for data_file in data_files:
            target_ligand = data_file.split("_G_CORR_H_")[-1]
            print("Loading " + data_file)
            cs[label][ref_ligand][target_ligand] = np.loadtxt(data_file)[0]
