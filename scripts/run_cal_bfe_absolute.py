"""
run calculations of binding free energies using the absolute (or apo) estimator
"""
from __future__ import print_function

import os
import argparse
import glob
import numpy as np
import copy
np.seterr(over='raise')     # raise exception if overload

from _absolute_estimators import MultiStruScores
from _algdock import load_algdock_snapshots_for_each_of_six_yank_systems


def averaging(yank_systems, result_dir, weights):
    """
    :param yank_systems: list of str
    :param result_dir: str
    :param weights: dict,
                    weights[ref_ligand_name][snapshot] -> float
                    weights["systems"][ref_ligand_name] -> float
    :return: None
    """
    #
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
            scores = MultiStruScores(args.scores_dir, group, code,  weights, yank_systems)
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
                    raise RuntimeError("unknown combining rule")

                id = scores.get_id()
                for FF in averages.keys():
                    out_files[rule][FF].write("%s   %20.10f %20.10f\n" %(id, averages[FF], standard_dev[FF]))

    for rule in combining_rules:
        for FF in FFs:
            out_files[rule][FF].close()

    return None


def equalize_system_weights(original_weights):
    """
    make all the ref systems equally weighted
    :param original_weights: dict
    :return: new_weights, dict
    """
    new_weights = copy.deepcopy(original_weights)
    for system in new_weights["systems"].keys():
        new_weights["systems"][system] = 1.0
    return new_weights


def take_6_apo(original_weights):
    """
    only snapshots from apo of YANK states get nonzero weights
    :param original_weights: dict
    :return: new_weights
    """
    new_weights = equalize_system_weights(original_weights)
    ordered_snapshots = load_algdock_snapshots_for_each_of_six_yank_systems() 

    for system in new_weights["systems"].keys():
        for i, snapshot in enumerate(ordered_snapshots[system]):
            if i < len(ordered_snapshots[system]) - 6:
                new_weights[system][snapshot] = 0.
            else:
                new_weights[system][snapshot] = 1.
            
    return new_weights


def take_12_near_apo(original_weights):
    """
    :param original_weights: dict
    :return: dict
    """
    new_weights = equalize_system_weights(original_weights)
    ordered_snapshots = load_algdock_snapshots_for_each_of_six_yank_systems()

    for system in new_weights["systems"].keys():
        for i, snapshot in enumerate( ordered_snapshots[system] ):
            if i < len(ordered_snapshots[system]) - 12:
                new_weights[system][snapshot] = 0.
    return new_weights


def take_24_near_apo(original_weights):
    """
    :param original_weights: dict
    :return: dict
    """
    new_weights = equalize_system_weights(original_weights)
    ordered_snapshots = load_algdock_snapshots_for_each_of_six_yank_systems()

    for system in new_weights["systems"].keys():
        for i, snapshot in enumerate(ordered_snapshots[system]):
            if i < len(ordered_snapshots[system]) - 24:
                new_weights[system][snapshot] = 0.
    return new_weights


DOCK6_SUB_DIR = "dock6"

parser = argparse.ArgumentParser()
parser.add_argument("--scores_dir", type=str,
    default="~/T4_Lysozyme/Bing_Calculations/Correct_Pro_BornRadii/Concatened_Scores/version2016_May_AlGDock/OBC2")


# "all96", "6apo", "12nearapo", "24nearapo"
parser.add_argument("--which_snapshots_to_take", type=str, default="all96")
args = parser.parse_args()

assert args.which_snapshots_to_take in ["all96", "6apo", "12nearapo", "24nearapo"], "unknown which_snapshots_to_take"
print("use " + args.which_snapshots_to_take)

combining_rules = ["Mean", "ExpMean", "Min"]

ligand_groups = glob.glob(os.path.join(args.scores_dir, DOCK6_SUB_DIR, "*"))
ligand_groups = [os.path.basename(d) for d in ligand_groups]

if os.path.basename(args.scores_dir) == "OBC2":
    print("OBC2 weights")
    from load_mbar_weights_apo_OBC2 import load_mbar_weights
elif os.path.basename(args.scores_dir) == "PBSA":
    print("PBSA weights")
    from load_mbar_weights_apo_PBSA import load_mbar_weights
else:
    raise ValueError("unknown phase "+os.path.basename(args.scores_dir))

ligand_3l_codes = {}
for group in ligand_groups:
    codes = glob.glob(os.path.join(args.scores_dir, DOCK6_SUB_DIR, group, "*"))
    ligand_3l_codes[group] = [os.path.basename(c) for c in codes]


block_weights, state_weights, single_snap_weights, stru_group_weights_equal_sys, stru_group_weights_ub_weighted = load_mbar_weights()

yank_systems = [key for key in block_weights.keys() if key not in ["systems"]]
print("yank systems ", yank_systems)


state_weights_equal_systems = equalize_system_weights(state_weights)

print(ligand_groups[0], ligand_3l_codes[ligand_groups[0]][0])
FFs = MultiStruScores(args.scores_dir, ligand_groups[0], ligand_3l_codes[ligand_groups[0]][0],
                      block_weights, yank_systems).get_FFs()
print(FFs)

if args.which_snapshots_to_take == "all96":
    use_state_weights = equalize_system_weights(state_weights)

elif args.which_snapshots_to_take == "6apo":
    use_state_weights = take_6_apo(state_weights)

elif args.which_snapshots_to_take == "12nearapo":
    use_state_weights = take_12_near_apo(state_weights)

elif args.which_snapshots_to_take == "24nearapo":
    use_state_weights = take_24_near_apo(state_weights)
else:
    raise ValueError("Unknown which_snapshots_to_take")

averaging(yank_systems, "all_snap__equal_sys__state_weight", use_state_weights)


