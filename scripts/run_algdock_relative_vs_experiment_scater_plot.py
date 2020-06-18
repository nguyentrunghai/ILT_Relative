"""
"""

from __future__ import print_function

import os
import argparse

import numpy as np

from _algdock import SIX_YANK_SYSTEMS
from _yank import load_scores, load_exper_bfes
from _plots import scatter_plot, scatter_plot_info

parser = argparse.ArgumentParser()

parser.add_argument("--exper_results", type=str, default="/home/tnguye46/T4_Lysozyme/Experiments/experiment_dG.dat")
parser.add_argument("--yank_results", type=str, default="/home/tnguye46/T4_Lysozyme/Yank/yank_results.dat")

# /home/tnguye46/T4_Lysozyme/Relative_Binding_FE/Relative_FE_Est_1/all96
# /home/tnguye46/T4_Lysozyme/Relative_Binding_FE/Relative_FE_Est_with_CV_method_4a/flip_sign_c__not_subtract_self
parser.add_argument("--algdock_dir", type=str, default="AlGDock")

# "with" or "without"
parser.add_argument("--control_variate", type=str, default=" ")

parser.add_argument("--error_scale_factor", type=float, default=1.)

parser.add_argument("--result_dir_suffix", type=str, default="__equal_sys__single_weight")
parser.add_argument("--combining_rule", type=str, default="ExpMean")
parser.add_argument("--rel_bfe_file", type=str, default="OpenMM_OBC2_MBAR.score")

parser.add_argument("--xlabel", type=str, default="Experimental results (kcal/mol)")
parser.add_argument("--ylabel", type=str, default="AlGDock free energy (kcal/mol)")

args = parser.parse_args()

assert args.control_variate in ["with", "without"], "Unknown control_variate"

# load experiment results
exper_bfes = load_exper_bfes(args.exper_results, id_col=1, score_col=2, exclude_ligands=[])

ref_ligands = SIX_YANK_SYSTEMS
print("Ref ligands:\n" + "\n".join(ref_ligands))
target_ligands = [lig for lig in exper_bfes.keys() if lig not in ref_ligands]
print("Target ligands:\n" + "\n".join(target_ligands))

# load yank results
yank_bfes, _ = load_scores(args.yank_results, 0, 1, 2, exclude_ligands=[])

# rbfe WITH CV
if args.control_variate == "with":
    algdock_bfes = {}            # algdock_bfes[ref_ligand][target_ligand] -> float
    algdock_errors = {}
    for ref_ligand in ref_ligands:
        infile = os.path.join(args.algdock_dir, ref_ligand + args.result_dir_suffix,
                              args.combining_rule, args.rel_bfe_file)
        print("Loading", infile)
        bfes, errors = load_scores(infile, 0, 1, 2, exclude_ligands=ref_ligands)

        # keep only target ligands
        bfes = {ligand: bfes[ligand] for ligand in bfes if ligand in target_ligands}
        errors = {ligand: errors[ligand] for ligand in errors if ligand in target_ligands}

        # add reference value from YANK
        for target_ligand in bfes:
            bfes[target_ligand] = bfes[target_ligand] + yank_bfes[ref_ligand]

        algdock_bfes[ref_ligand] = bfes
        algdock_errors[ref_ligand] = errors

# rbfe WITHOUT CV
if args.control_variate == "without":
    short_names = ["active.D__DAA", "active.C__CAA", "active.C__CAB"]
    full_names = ["lysozyme." + lig for lig in short_names]

    algdock_bfes = {}            # algdock_bfes[ref_ligand][target_ligand] -> float
    algdock_errors = {}
    for ref_ligand in ref_ligands:
        infile = os.path.join(args.algdock_dir, ref_ligand + args.result_dir_suffix,
                              args.combining_rule, args.rel_bfe_file)
        print("Loading", infile)
        bfes, errors = load_scores(infile, 0, 1, 2, exclude_ligands=[])

        # change name for some ligands
        for s_name, f_name in zip(short_names, full_names):
            bfes[f_name] = bfes[s_name]
            errors[f_name] = errors[s_name]

        self_bfe = bfes[ref_ligand]

        # keep only target ligands
        bfes = {ligand: bfes[ligand] for ligand in bfes if ligand in target_ligands}
        errors = {ligand: errors[ligand] for ligand in errors if ligand in target_ligands}

        for target_ligand in bfes:
            bfes[target_ligand] = bfes[target_ligand] - self_bfe + yank_bfes[ref_ligand]

        algdock_bfes[ref_ligand] = bfes
        algdock_errors[ref_ligand] = errors

# plot
xs = []
for target_ligand in target_ligands:
    xs.append(exper_bfes[target_ligand])
xs = np.array(xs)

for ref_ligand in ref_ligands:
    out_dir = ref_ligand
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_fig = os.path.join(out_dir, "algdock_vs_experiment.pdf")
    out_log = os.path.join(out_dir, "rmse_pearsonR.dat")

    ys = []
    y_errs = []
    for target_ligand in target_ligands:
        ys.append(algdock_bfes[ref_ligand][target_ligand])
        y_errs.append(algdock_errors[ref_ligand][target_ligand])

    ys = np.array(ys)
    y_errs = np.array(y_errs) / 2.

    dummy_ligands = ["abc" for _ in ys]
    scatter_plot_info(xs, ys, dummy_ligands, out_log)
    scatter_plot(xs, ys, args.xlabel, args.ylabel, "without_CV.pdf",
                 show_xy_axes=True,
                 yerr=y_errs,
                 show_regression_line=True,
                 show_diagonal_line=False,
                 show_rmse=True,
                 show_R=True,
                 show_regression_line_eq=True,
                 markersize=4,
                 same_xy_scale=False,
                 text_pos=[0.1, 0.7])


print("DONE")
