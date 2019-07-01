"""
"""

from __future__ import print_function

import os
import argparse

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

from _yank import load_scores
from _algdock import SIX_YANK_SYSTEMS

parser = argparse.ArgumentParser()

parser.add_argument("--yank_results", type=str, default="/home/tnguye46/T4_Lysozyme/Yank/yank_results.dat")

parser.add_argument("--rel_bfe_dir", type=str,
default="/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/Relative_FE_Est_with_CV_using_exp_mean_method_3/with_neg_cap/all96")
parser.add_argument("--result_dir_suffix", type=str, default="__equal_sys__single_weight")
parser.add_argument("--combining_rule", type=str, default="ExpMean")
parser.add_argument("--rel_bfe_file", type=str, default="OpenMM_OBC2_MBAR.score")

args = parser.parse_args()

ref_ligands = SIX_YANK_SYSTEMS
yank_bfes, yank_bfe_errors = load_scores(args.yank_results, 0, 1, 2, [])

all_ligands = yank_bfes.keys()
target_ligands = [ligand for ligand in all_ligands if ligand not in ref_ligands]

# deviation from yank
bfe_devs = {}            # bfe_devs[ref_ligand][target_ligand] -> float
bfe_errors = {}
for ref_ligand in ref_ligands:

    infile = os.path.join(args.rel_bfe_dir, ref_ligand + args.result_dir_suffix,
                          args.combining_rule, args.rel_bfe_file)
    print("Loading", infile)
    bfes, errors = load_scores(infile, 0, 1, 2, ref_ligands)

    for target_ligand in bfes:
        bfes[target_ligand] = bfes[target_ligand] + yank_bfes[ref_ligand] - yank_bfes[target_ligand]

    bfe_devs[ref_ligand] = bfes
    bfe_errors[ref_ligand] = errors

# load C and g-h correlation coefficients
cs = {}
c_errors = {}
corr_coeffs = {}
corr_errors = {}

for ref_ligand in ref_ligands:
    cs[ref_ligand] = {}
    c_errors[ref_ligand] = {}
    corr_coeffs[ref_ligand] = {}
    corr_errors[ref_ligand] = {}

    for target_ligand in target_ligands:
        infile = os.path.join(args.rel_bfe_dir, ref_ligand + args.result_dir_suffix,
                              args.combining_rule, ref_ligand + "_G_CORR_H_" + target_ligand)
        print("Loading", infile)
        c, c_err, corr, corr_err = np.loadtxt(infile)

        cs[ref_ligand][target_ligand] = c
        c_errors[ref_ligand][target_ligand] = c_err
        corr_coeffs[ref_ligand][target_ligand] = corr
        corr_errors[ref_ligand][target_ligand] = corr_err


FONTSIZE = 8
FONT = {"fontname": "Arial"}

# plot bfe_devs vs cs
xs = [cs[ref_ligand][target_ligand] for ref_ligand in ref_ligands for target_ligand in target_ligands]
ys = [bfe_devs[ref_ligand][target_ligand] for ref_ligand in ref_ligands for target_ligand in target_ligands]

xs = np.array(xs)
ys = np.array(ys)

fig, ax = plt.subplots(1, 1, figsize=(3.2, 2.4))
ax.scatter(xs, ys)
# scilimits: (m, n), pair of integers; if style is 'sci', scientific notation will be used for
# numbers outside the range 10**m to 10**n. Use (0,0) to include all numbers.

ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
ax.set_xlabel("$C$", fontsize=FONTSIZE, **FONT)
ax.set_ylabel("Deviation from YANK (kcal/mol)", fontsize=FONTSIZE, **FONT)

fig.tight_layout()
fig.savefig("bfe_dev_vs_C.pdf")

# remove big outliers
idx = xs.argsort()
xs_1 = xs[idx][1:-15]
ys_1 = ys[idx][1:-15]

fig, ax = plt.subplots(1, 1, figsize=(3.2, 2.4))
ax.scatter(xs_1, ys_1)

ax.set_xlabel("$C$", fontsize=FONTSIZE, **FONT)
ax.set_ylabel("Deviation from YANK (kcal/mol)", fontsize=FONTSIZE, **FONT)

fig.tight_layout()
fig.savefig("bfe_dev_vs_C_remove_big_values.pdf")

#-----------------------------------------------
# plot bfe_devs vs c_errors
xs = [c_errors[ref_ligand][target_ligand] for ref_ligand in ref_ligands for target_ligand in target_ligands]
xs = np.array(xs)
xs_1 = xs[idx][1:-15]

fig, ax = plt.subplots(1, 1, figsize=(3.2, 2.4))
ax.scatter(xs_1, ys_1)

ax.set_xlabel("Error in $C$", fontsize=FONTSIZE, **FONT)
ax.set_ylabel("Deviation from YANK (kcal/mol)", fontsize=FONTSIZE, **FONT)

fig.tight_layout()
fig.savefig("bfe_dev_vs_C_errors.pdf")

#--------------------------------------------------
# plot bfe_devs vs corr
xs = [corr_coeffs[ref_ligand][target_ligand] for ref_ligand in ref_ligands for target_ligand in target_ligands]
xs = np.array(xs)

fig, ax = plt.subplots(1, 1, figsize=(3.2, 2.4))
ax.scatter(xs, ys)

ax.set_xlabel("Corr($g, h$)", fontsize=FONTSIZE, **FONT)
ax.set_ylabel("Deviation from YANK (kcal/mol)", fontsize=FONTSIZE, **FONT)

fig.tight_layout()
fig.savefig("bfe_dev_vs_corr.pdf")

#---------------------------------------------------
# plot bfe_devs vs corr_error
xs = [corr_errors[ref_ligand][target_ligand] for ref_ligand in ref_ligands for target_ligand in target_ligands]
xs = np.array(xs)

fig, ax = plt.subplots(1, 1, figsize=(3.2, 2.4))
ax.scatter(xs, ys)

ax.set_xlabel("Error in Corr($g, h$)", fontsize=FONTSIZE, **FONT)
ax.set_ylabel("Deviation from YANK (kcal/mol)", fontsize=FONTSIZE, **FONT)

fig.tight_layout()
fig.savefig("bfe_dev_vs_corr_error.pdf")