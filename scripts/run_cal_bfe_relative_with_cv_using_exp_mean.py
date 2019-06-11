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
from _relative_estimators import relative_bfe_with_cv_using_exp_mean

parser = argparse.ArgumentParser()

parser.add_argument("--scores_dir", type=str,
                    default="/home/tnguye46/T4_Lysozyme/Bing_Calculations/Correct_Pro_BornRadii/Concatened_Scores/version2016_May_AlGDock/OBC2")

parser.add_argument("--interaction_energies_dir", type=str,
                    default="/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/OpenMM_OBC2_interaction_energies_for_576_algdock_snapshots")

parser.add_argument("--FF", type=str, default="OpenMM_OBC2_MBAR")

parser.add_argument("--bootstrap_repeats", type=int, default=1000)

parser.add_argument("--result_dir_suffix", type=str, default="__equal_sys__single_weight")

parser.add_argument("--combining_rule", type=str, default="ExpMean")

args = parser.parse_args()

_, _, single_snap_weights, _, _ = load_mbar_weights()
ref_ligands = [ligand for ligand in single_snap_weights.keys() if ligand != "systems"]
print("ref_ligands", ref_ligands)

yank_interaction_energies = load_interaction_energies(path=args.interaction_energies_dir)

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
        hs, gs, rel_bfe = relative_bfe_with_cv_using_exp_mean(snapshots, args.scores_dir, target_ligand, ref_ligand,
                                            single_snap_weights, yank_interaction_energies, args.FF, verbose=True)

        bootstrap_ests = []
        for _ in range(args.bootstrap_repeats):
            random_snapshots = np.random.choice(snapshots, size=len(snapshots), replace=True)
            _, _, bfe = relative_bfe_with_cv_using_exp_mean(random_snapshots, args.scores_dir, target_ligand,
                                                            ref_ligand, single_snap_weights, yank_interaction_energies,
                                                            args.FF, verbose=False)
            if (not np.isnan(bfe)) and (not np.isinf(bfe)):
                bootstrap_ests.append(bfe)

        error = np.std(bootstrap_ests)

        out_file_handle.write("%s   %20.10f %20.10f\n" %(target_ligand, rel_bfe, error))

        rel_bfe_vs_self = os.path.join(result_dir, args.combining_rule, ref_ligand + "_vs_" + target_ligand)

        with open(rel_bfe_vs_self, "w") as handle:
            handle.write("# h          g\n")
            for h, g in zip(hs, gs):
                handle.write("%20.10f %20.10f\n" % (h, g))

    out_file_handle.close()

