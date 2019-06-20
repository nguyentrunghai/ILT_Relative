"""
"""

import os
import argparse

from _yank import load_scores
from _algdock import SIX_YANK_SYSTEMS

parser = argparse.ArgumentParser()

parser.add_argument("--yank_results", type=str, default="/home/tnguye46/T4_Lysozyme/Yank/yank_results.dat")

parser.add_argument("--rel_bfe_dir", type=str,
default="/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/Relative_FE_Est_with_CV_using_exp_mean_with_neg_cap/all96")
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




