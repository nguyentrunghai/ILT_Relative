"""
Run calculations of relative binding free energies using holo estimator
"""
from __future__ import print_function

import os
import argparse
import glob
import numpy as np

from _process_yank_outputs import load_interaction_energies

from _relative_estimators import RelBFEWithoutCV

from _weight_processing import equalize_system_weights, take_6_holo, take_12_near_holo, take_24_near_holo

np.seterr(over='raise')     # raise exception if overload


def averaging(scores_dir, ligand_3l_codes, yank_systems, result_dir,
              weights, yank_interaction_energies,
              combining_rules):
    """
    :param scores_dir: str
    :param ligand_3l_codes: dict, ligand_3l_codes[group] -> list of three-letter strings
    :param yank_systems: list of str
    :param result_dir: str
    :param weights: dict,
                    weights[ref_ligand_name][snapshot] -> float
                    weights["systems"][ref_ligand_name] -> float
    :param yank_interaction_energies: dict, yank_interaction_energies[system][snapshot] -> float
    :param combining_rules: list of str
    :return: None
    """
    if not os.path.isdir(result_dir):
        os.system("mkdir " + result_dir)

    out_files = {}
    for rule in combining_rules:
        if not os.path.isdir(os.path.join(result_dir, rule)):
            os.system("mkdir " + os.path.join(result_dir, rule))
        out_files[rule] = {}
        for FF in FFs:
            out_files[rule][FF] = open(os.path.join(result_dir, rule, FF + ".score"), "w")

    for group in ligand_3l_codes.keys():

        for code in ligand_3l_codes[group]:
            scores = RelBFEWithoutCV(scores_dir, group, code,  weights, yank_systems, yank_interaction_energies)
            scores.check_extreme_low()

            for rule in combining_rules:
                if rule == "Mean":
                    averages = scores.get_mean()
                    standard_dev = scores.get_mean_std()

                elif rule == "Min":
                    averages = scores.get_min()
                    standard_dev = scores.get_min_std()

                elif rule == "ExpMean":
                    averages = scores.get_exp_mean()
                    standard_dev = scores.get_exp_mean_std()
                else:
                    raise ValueError("unknown combining rule")

                id = scores.get_id()
                for FF in averages.keys():
                    out_files[rule][FF].write("%s   %20.10f %20.10f\n" %(id, averages[FF], standard_dev[FF]))

    for rule in combining_rules:
        for FF in FFs:
            out_files[rule][FF].close()

    return None


DOCK6_SUB_DIR = "dock6"

parser = argparse.ArgumentParser()

parser.add_argument("--scores_dir", type=str,
default="~/T4_Lysozyme/Bing_Calculations/Correct_Pro_BornRadii/Concatened_Scores/version2016_May_AlGDock/OBC2")

parser.add_argument("--interaction_energies_dir", type=str,
default="~/T4_Lysozyme/Relative_Binding_FE/OpenMM_OBC2_interaction_energies_for_576_algdock_snapshots")

# "all96", "6holo", "12nearholo", "24nearholo"
parser.add_argument("--which_snapshots_to_take", type=str, default="all96")

args = parser.parse_args()

assert args.which_snapshots_to_take in ["all96", "6holo", "12nearholo", "24nearholo"], "unknown which_snapshots_to_take"
print("use", args.which_snapshots_to_take)

combining_rules = ["Mean", "ExpMean", "Min"]

ligand_groups = glob.glob(os.path.join(args.scores_dir, DOCK6_SUB_DIR, "*"))
ligand_groups = [os.path.basename(d) for d in ligand_groups]

if os.path.basename(args.scores_dir) == "OBC2":
    print("OBC2 weights")
    from load_mbar_weights_holo_OBC2 import load_mbar_weights
elif os.path.basename(args.scores_dir) == "PBSA":
    print("PBSA weights")
    from load_mbar_weights_holo_PBSA import load_mbar_weights
else:
    raise ValueError("unknown phase "+os.path.basename(args.scores_dir))


ligand_3l_codes = {}
for group in ligand_groups:
    codes = glob.glob(os.path.join( args.scores_dir, DOCK6_SUB_DIR, group, "*"))
    ligand_3l_codes[group] = [os.path.basename(c) for c in codes]


block_weights, state_weights, single_snap_weights, stru_group_weights_equal_sys, stru_group_weights_ub_weighted = load_mbar_weights()

yank_interaction_energies = load_interaction_energies(path=args.interaction_energies_dir) 

yank_systems = [key for key in block_weights.keys() if key not in ["systems"]]
print("yank systems ", yank_systems)

state_weights_equal_systems = equalize_system_weights(state_weights)
single_snap_weights_equal_systems = equalize_system_weights(single_snap_weights) 

print(ligand_groups[0], ligand_3l_codes[ligand_groups[0]][0])

FFs = RelBFEWithoutCV(args.scores_dir, ligand_groups[0], ligand_3l_codes[ligand_groups[0]][0],
                      block_weights, yank_systems, yank_interaction_energies).get_FFs()
print(FFs)

if args.which_snapshots_to_take == "all96":

    use_state_weights = equalize_system_weights(state_weights)
    use_single_weights = equalize_system_weights(single_snap_weights)

elif args.which_snapshots_to_take == "6holo":
    use_state_weights = take_6_holo(state_weights)
    use_single_weights = take_6_holo(single_snap_weights)

elif args.which_snapshots_to_take == "12nearholo":
    use_state_weights = take_12_near_holo(state_weights)
    use_single_weights = take_12_near_holo(single_snap_weights)

elif args.which_snapshots_to_take == "24nearholo":
    use_state_weights = take_24_near_holo(state_weights)
    use_single_weights = take_24_near_holo(single_snap_weights)
else:
    raise ValueError("")

for y_sys in yank_systems:
    # a yank system is a ref ligand
    averaging(args.scores_dir, ligand_3l_codes, [y_sys],
              y_sys + "__equal_sys__state_weight",
              use_state_weights, yank_interaction_energies, combining_rules)


    averaging(args.scores_dir, ligand_3l_codes, [y_sys],
              y_sys + "__equal_sys__single_weight",
              use_single_weights, yank_interaction_energies, combining_rules)
