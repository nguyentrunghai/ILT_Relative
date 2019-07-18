"""
plot histograms of c for individual methods and scatter plots for pairs of methods
"""

from __future__ import print_function

import argparse
import glob
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

from _algdock import SIX_YANK_SYSTEMS

parser = argparse.ArgumentParser()
parser.add_argument("--data_dirs", type=str,
                    default="../Relative_FE_Est_with_CV_method_2a/not_flip_sign_c ../Relative_FE_Est_with_CV_method_2b/not_flip_sign_c")
parser.add_argument("--method_labels", type=str,
                    default="method_2a method_2b")
parser.add_argument("--threshold_4_outliers", type=float, default=5)
args = parser.parse_args()


def _is_outlier(arr, test_value, threshold):
    assert threshold > 0, "threshold must be positive"
    iqr = np.percentile(arr, 75.) - np.percentile(arr, 25.)
    score = np.abs(test_value - np.median(arr))
    return score > threshold*iqr


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
        data_files = [data_file for data_file in data_files if data_file.split("_G_CORR_H_")[-1] not in ref_ligands]

        for data_file in data_files:
            target_ligand = data_file.split("_G_CORR_H_")[-1]
            print("Loading " + data_file)
            cs[label][ref_ligand][target_ligand] = np.loadtxt(data_file)[0]

        target_ligands = cs[label][ref_ligand].keys()

for label_x in method_labels:
    for label_y in method_labels:
        if label_y != label_x:
            xs = []
            ys = []

            for ref_ligand in ref_ligands:
                for target_ligand in target_ligands:
                    xs.append(cs[label_x][ref_ligand][target_ligand])
                    ys.append(cs[label_y][ref_ligand][target_ligand])

            xs_new = []
            ys_new = []
            for x, y in zip(xs, ys):
                if not _is_outlier(xs, x, args.threshold_4_outliers):
                    if not _is_outlier(ys, y, args.threshold_4_outliers):
                        xs_new.append(x)
                        ys_new.append(y)

            xs_new = np.array(xs_new)
            ys_new = np.array(ys_new)

            data = pd.DataFrame({label_x: xs_new, label_y: ys_new})
            #fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3.2, 2.4))
            plt.figure(figsize=(3.2, 2.4))
            sns.jointplot(x=label_x, y=label_y, data=data, kind="scatter")
            plt.xlabel(label_x)
            plt.ylabel(label_y)
            plt.tight_layout()

            out = label_x + "_vs_" + label_y + ".pdf"
            dpi = 300
            plt.savefig(out, dpi=dpi)
