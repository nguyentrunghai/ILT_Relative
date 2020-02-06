"""
"""

from __future__ import print_function

import os
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

from _yank import load_scores
from _algdock import SIX_YANK_SYSTEMS
from _plots import scatter_plot, scatter_plot_info

parser = argparse.ArgumentParser()

parser.add_argument("--yank_results", type=str, default="/home/tnguye46/T4_Lysozyme/Yank/yank_results.dat")

parser.add_argument("--bfe_without_cv_dir", type=str, default="Relative_FE_Est_1/all96")

parser.add_argument("--bfe_with_cv_dir", type=str, default="flip_sign_c__not_subtract_self")

parser.add_argument("--error_scale_factor", type=float, default=1.)

parser.add_argument("--result_dir_suffix", type=str, default="__equal_sys__single_weight")
parser.add_argument("--combining_rule", type=str, default="ExpMean")
parser.add_argument("--rel_bfe_file", type=str, default="OpenMM_OBC2_MBAR.score")

args = parser.parse_args()

# load yank results
ref_ligands = SIX_YANK_SYSTEMS
yank_bfes, yank_bfe_errors = load_scores(args.yank_results, 0, 1, 2, exclude_ligands=[])

all_ligands = yank_bfes.keys()
target_ligands = [ligand for ligand in all_ligands if ligand not in ref_ligands]


# rbfe WITH CV
rbfes_with_cv = {}            # rbfes_with_cv[ref_ligand][target_ligand] -> float
rbfe_errors_with_cv = {}
for ref_ligand in ref_ligands:
    infile = os.path.join(args.bfe_with_cv_dir, ref_ligand + args.result_dir_suffix,
                          args.combining_rule, args.rel_bfe_file)
    print("Loading", infile)
    bfes, errors = load_scores(infile, 0, 1, 2, exclude_ligands=ref_ligands)

    for target_ligand in bfes:
        bfes[target_ligand] = bfes[target_ligand] + yank_bfes[ref_ligand]

    rbfes_with_cv[ref_ligand] = bfes
    rbfe_errors_with_cv[ref_ligand] = errors


# rbfe WITHOUT CV
rbfes_without_cv = {}            # rbfes_without_cv[ref_ligand][target_ligand] -> float
rbfe_errors_without_cv = {}
for ref_ligand in ref_ligands:
    infile = os.path.join(args.bfe_without_cv_dir, ref_ligand + args.result_dir_suffix,
                          args.combining_rule, args.rel_bfe_file)
    print("Loading", infile)
    bfes, errors = load_scores(infile, 0, 1, 2, exclude_ligands=[])
    self_bfe = bfes[ref_ligand]

    # remove reference ligands in the dict
    bfes = {ligand: bfes[ligand] for ligand in bfes if ligand in target_ligands}
    errors = {ligand: errors[ligand] for ligand in errors if ligand in target_ligands}

    for target_ligand in bfes:
        bfes[target_ligand] = bfes[target_ligand] - self_bfe + yank_bfes[ref_ligand]

    rbfes_without_cv[ref_ligand] = bfes
    rbfe_errors_without_cv[ref_ligand] = errors

# plot
xs = []
x_errs = []
for ref_ligand in ref_ligands:
    for target_ligand in target_ligands:
        xs.append(yank_bfes[target_ligand])
        x_errs.append(yank_bfe_errors[target_ligand])

xs = np.array(xs)
x_errs = np.aray(x_errs) / 2.

# without CV
ys = []
y_errs = []
for ref_ligand in ref_ligands:
    for target_ligand in target_ligands:
        ys.append(rbfes_without_cv[ref_ligand][target_ligand])
        y_errs.append(rbfe_errors_without_cv[ref_ligand][target_ligand])

ys = np.array(ys)
y_errs = np.aray(y_errs) / 2.



FONTSIZE = 8
FONT = {"fontname": "Arial"}

# plot difference in absolute deviation vs corr coef
xs = [corr_coeffs[ref_ligand][target_ligand] for ref_ligand in ref_ligands for target_ligand in target_ligands]
ys = [np.abs(devs_with_cv[ref_ligand][target_ligand]) - np.abs(devs_without_cv[ref_ligand][target_ligand])
      for ref_ligand in ref_ligands for target_ligand in target_ligands]

xs = np.array(xs)
ys = np.array(ys) - args.shift

fig, ax = plt.subplots(1, 1, figsize=(3.2, 2.4))
ax.scatter(xs, ys)
ax.axhline(y=0, c="k")

ax.set_xlabel("Corr($g, h$)", fontsize=FONTSIZE, **FONT)
ax.set_ylabel("Diff. in Abs. Dev. from YANK (kcal/mol)", fontsize=FONTSIZE, **FONT)

fig.tight_layout()
fig.savefig("dev_diff_vs_corr.pdf")


# rate of negative Diff. in Abs. Dev. from YANK grouped by correlation bin
rate_neg_diff_dev = bin_corr_coef(xs, ys)
rate_neg_diff_dev.to_csv("rate_neg_diff_dev.csv")

overall_rate_neg_diff_dev = (ys < 0.).mean()
print("Overall rate of negative difference in absolute YANK deviation: %0.5f" % overall_rate_neg_diff_dev)


# plot estimation errors with vs without cv
xs = [bfe_errors_without_cv[ref_ligand][target_ligand] for ref_ligand in ref_ligands
      for target_ligand in target_ligands]
ys = [bfe_errors_with_cv[ref_ligand][target_ligand] for ref_ligand in ref_ligands
      for target_ligand in target_ligands]

xs = np.array(xs)
ys = np.array(ys) * args.error_scale_factor

fig, ax = plt.subplots(1, 1, figsize=(3.2, 2.4))
ax.scatter(xs, ys)

lower = np.min([xs.min(), xs.min()])
upper = np.max([xs.max(), xs.max()])

ax.plot([lower, upper], [lower, upper], c="k")
ax.set_xlim([lower, upper])
ax.set_ylim([lower, upper])

ax.set_xlabel("Bootstrap std without CV (kcal/mol)", fontsize=FONTSIZE, **FONT)
ax.set_ylabel("Bootstrap std with CV (kcal/mol)", fontsize=FONTSIZE, **FONT)

fig.tight_layout()
fig.savefig("error_with_vs_without_CV.pdf")

rate_error_with_less_than_without = (ys < xs).mean()
print("Rate at which errors of with are less than without CV %0.5f" % rate_error_with_less_than_without)

print("Mean bootstrap std without CV %0.5f" % xs.mean())
print("Mean bootstrap std with CV %0.5f" % ys.mean())
