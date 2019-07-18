"""
plot histograms of c for individual methods and scatter plots for pairs of methods
"""

from __future__ import print_function

import argparse
import glob
import os

from _algdock import SIX_YANK_SYSTEMS

parser = argparse.ArgumentParser()
parser.add_argument("--data_dirs", type=str,
                    default="../Relative_FE_Est_with_CV_method_2a/not_flip_sign_c")
parser.add_argument("--method_labels", type=str,
                    default="method_2a")
args = parser.parse_args()

data_dirs = args.data_dirs.split()
method_labels = args.method_labels.split()

ref_ligands = SIX_YANK_SYSTEMS
SUB_DIR_SUFFIX = "__equal_sys__single_weight"

cs = {}
for label, data_dir in zip(method_labels, data_dirs):
    cs[label] = {}

    for ref_ligand in ref_ligands:
        data_files = glob.glob(os.path.join(data_dir, ref_ligand+SUB_DIR_SUFFIX, "ExpMean", ref_ligand+"_G_CORR_H_*"))
        data_files = [f for f in data_files if f.split("_G_CORR_H_")[-1] != ref_ligand]
        print(data_files)