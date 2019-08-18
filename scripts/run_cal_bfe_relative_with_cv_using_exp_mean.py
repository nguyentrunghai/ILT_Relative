"""
run calculations for relative binding free energies using the holo estimator with
control variates
"""
from __future__ import print_function

import os
import argparse

import numpy as np

from _yank import YANK_LIGANDS
from load_mbar_weights_holo_OBC2 import load_mbar_weights
from _process_yank_outputs import load_interaction_energies
from _relative_estimators import relative_bfe_with_cv_using_exp_mean_method_2a
from _relative_estimators import relative_bfe_with_cv_using_exp_mean_method_2b
from _relative_estimators import relative_bfe_with_cv_using_exp_mean_method_3a
from _relative_estimators import relative_bfe_with_cv_using_exp_mean_method_3b

parser = argparse.ArgumentParser()

parser.add_argument("--scores_dir", type=str,
                    default="/home/tnguye46/T4_Lysozyme/Bing_Calculations/Correct_Pro_BornRadii/Concatened_Scores/version2016_May_AlGDock/OBC2")

parser.add_argument("--interaction_energies_dir", type=str,
                    default="/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/OpenMM_OBC2_interaction_energies_for_576_algdock_snapshots")

parser.add_argument("--FF", type=str, default="OpenMM_OBC2_MBAR")

parser.add_argument("--bootstrap_repeats", type=int, default=1000)

parser.add_argument("--result_dir_suffix", type=str, default="__equal_sys__single_weight")

parser.add_argument("--combining_rule", type=str, default="ExpMean")

# method is either 2a, 2b, 3a, 3b
parser.add_argument("--method", type=str, default="None")

parser.add_argument("--subtract_self", action="store_true", default=False)
parser.add_argument("--flip_sign_c", action="store_true", default=False)

args = parser.parse_args()


def _filter_big_values(arr):
    """
    Filter out big values so that np.sd does not give overflow error
    :param arr: 1d array or list
    :return: 1d array
    """
    arr = np.array(arr)
    arr_cen = arr - arr.mean()
    arr_out = []
    for v in arr_cen:
        try:
            v*v
        except FloatingPointError:
            print("Overflow if squared: %0.5e" % v)
        else:
            arr_out.append(v)
    return np.array(arr_out)


if args.method == "2a":
    print("Method 2a")
    relative_bfe_with_cv_using_exp_mean = relative_bfe_with_cv_using_exp_mean_method_2a

elif args.method == "2b":
    print("Method 2b")
    relative_bfe_with_cv_using_exp_mean = relative_bfe_with_cv_using_exp_mean_method_2b

elif args.method == "3a":
    print("Method 3a")
    relative_bfe_with_cv_using_exp_mean = relative_bfe_with_cv_using_exp_mean_method_3a

elif args.method == "3b":
    print("Method 3b")
    relative_bfe_with_cv_using_exp_mean = relative_bfe_with_cv_using_exp_mean_method_3b

else:
    raise ValueError("unrecognized method " + args.method)

_, _, single_snap_weights, _, _ = load_mbar_weights()
ref_ligands = [ligand for ligand in single_snap_weights.keys() if ligand != "systems"]
print("ref_ligands", ref_ligands)

yank_interaction_energies = load_interaction_energies(path=args.interaction_energies_dir)

if args.flip_sign_c:
    print("Flip sign of C if m_bar is negative")
else:
    print("Don't flip sign of C")

if args.subtract_self:
    print("Subtract from self relative BFE")

for ref_ligand in ref_ligands:
    print("Processing reference ligand", ref_ligand)

    result_dir = ref_ligand + args.result_dir_suffix
    if not os.path.isdir(result_dir):
        os.mkdir(result_dir)

    if not os.path.isdir(os.path.join(result_dir, args.combining_rule)):
        os.mkdir(os.path.join(result_dir, args.combining_rule))

    out_file = os.path.join(result_dir, args.combining_rule, args.FF + ".score")
    out_file_handle = open(out_file, "w")

    snapshots = single_snap_weights[ref_ligand].keys()

    for target_ligand in YANK_LIGANDS:
        hs, gs, c, correlation, rel_bfe = relative_bfe_with_cv_using_exp_mean(snapshots, args.scores_dir,
                                                                              target_ligand, ref_ligand,
                                                                              single_snap_weights,
                                                                              yank_interaction_energies, args.FF,
                                                                              subtract_self=args.subtract_self,
                                                                              flip_sign_c=args.flip_sign_c,
                                                                              verbose=True)

        bootstrap_bfes = []
        bootstrap_cs = []
        bootstrap_corrs = []
        for _ in range(args.bootstrap_repeats):
            random_snapshots = np.random.choice(snapshots, size=len(snapshots), replace=True)
            _, _, b_c, b_corr, bfe = relative_bfe_with_cv_using_exp_mean(random_snapshots, args.scores_dir,
                                                                         target_ligand, ref_ligand,
                                                                         single_snap_weights,
                                                                         yank_interaction_energies, args.FF,
                                                                         subtract_self=args.subtract_self,
                                                                         flip_sign_c=args.flip_sign_c,
                                                                         verbose=False)

            if (not np.isnan(bfe)) and (not np.isinf(bfe)):
                bootstrap_bfes.append(bfe)
                bootstrap_cs.append(b_c)
                bootstrap_corrs.append(b_corr)

        bootstrap_bfes = _filter_big_values(bootstrap_bfes)
        bootstrap_cs = _filter_big_values(bootstrap_cs)
        bootstrap_corrs = _filter_big_values(bootstrap_corrs)

        error = np.std(bootstrap_bfes)
        out_file_handle.write("%s   %20.10f %20.10f\n" %(target_ligand, rel_bfe, error))

        g_vs_h_out_file = os.path.join(result_dir, args.combining_rule, ref_ligand + "_G_VERSUS_H_" + target_ligand)
        with open(g_vs_h_out_file, "w") as handle:
            handle.write("# g          h\n")
            for g, h in zip(gs, hs):
                handle.write("%20.10e %20.10e\n" % (g, h))

        g_corr_h_out_file = os.path.join(result_dir, args.combining_rule, ref_ligand + "_G_CORR_H_" + target_ligand)
        with open(g_corr_h_out_file, "w") as handle:
            handle.write("# C        C_error      correlation    corr_error\n")
            handle.write("%20.10e %20.10e %20.10e %20.10e\n" % (c, np.std(bootstrap_cs),
                                                                correlation, np.std(bootstrap_corrs)))

    out_file_handle.close()

