"""
run convergence with respect to final for relative binding free energies without control variates
"""
from __future__ import print_function

import os
import argparse

import numpy as np

from _process_yank_outputs import load_interaction_energies
from load_mbar_weights_holo_OBC2 import load_mbar_weights
from _yank import load_scores
from _relative_estimators import RelBFEWithoutCV, relative_bfe_with_cv_using_exp_mean_method_3a
from _yank import YANK_LIGANDS

parser = argparse.ArgumentParser()
parser.add_argument("--algdock_score_dir", type=str,
                    default="Bing_Calculations/Correct_Pro_BornRadii/Concatened_Scores/version2016_May_AlGDock/OBC2")
parser.add_argument("--interaction_energies_dir",   type=str,
                    default="Relative_Binding_FE/OpenMM_OBC2_interaction_energies_for_576_algdock_snapshots")

parser.add_argument("--FF", type=str, default="OpenMM_OBC2_MBAR")
parser.add_argument("--weight_scheme", type=str, default="__equal_sys__single_weight")
parser.add_argument("--combining_rule", type=str, default="ExpMean")

parser.add_argument("--final_results_dir", type=str,
                    default="Relative_Binding_FE/Relative_FE_Est_1/all96")

parser.add_argument("--bootstrap_repeats", type=int, default=100)
parser.add_argument("--sample_sizes", type=str, default="10 50 96")

# "without_cv" or "with_cv_3a"
parser.add_argument("--method", type=str, default="none")

args = parser.parse_args()


def _load_final_fe(final_results_dir, ref_ligands, weight_scheme, combining_rule, FF):
    """
    :param final_results_dir: str
    :param ref_ligands: list of str
    :param weight_scheme: str
    :param combining_rule: str
    :param score_file: str
    :return final_fes: dict {ref_ligand (str): {ligand (str): fe (float)}}
    """
    final_fes = {}
    for ref_ligand in ref_ligands:
        fe_file = os.path.join(final_results_dir, ref_ligand+weight_scheme, combining_rule, FF+".score")
        print("Loading file for final fes:", fe_file)
        fes, _ = load_scores(fe_file, 0, 1, 2, exclude_ligands=[])
        fes.pop(ref_ligand)

        final_fes[ref_ligand] = fes
    return final_fes


def _pearson_r_rmse(reference_vals, target_vals):
    """
    :param reference_vals: dict, {ligand (str): free energy (float)}
    :param target_vals: dict, {ligand (str): free energy (float)}
    :return r: float
    """
    ligands = set(reference_vals.keys()).intersection(target_vals.keys())
    ligands = [ligand for ligand in ligands if str(reference_vals[ligand]).lower() not in ["inf", "-inf", "nan"]]
    ligands = [ligand for ligand in ligands if str(target_vals[ligand]).lower() not in ["inf", "-inf", "nan"]]

    xs = np.array([reference_vals[ligand] for ligand in ligands], dtype=float)
    ys = np.array([target_vals[ligand] for ligand in ligands], dtype=float)

    r = np.corrcoef([xs, ys])[0, -1]
    rmse = ((xs - ys) ** 2).mean()
    rmse = np.sqrt(rmse)

    return r, rmse


def _r_rmse_one_ref_ligand_a_random_sample_of_snapshot_without_cv(algdock_score_dir, target_ligands,
                                                      ref_ligand, ref_ligands,
                                                      FF, weights, yank_interaction_energies,
                                                      sample_size, final_fes):
    """
    :param algdock_score_dir: str
    :param target_ligands: list of str
    :param ref_ligand: str
    :param FF: str
    :param ref_ligands: list of str
    :param weights: dict,
                    weights[ref_ligand_name][snapshot] -> float
                    weights["systems"][ref_ligand_name] -> float
    :param yank_interaction_energies: dict, yank_interaction_energies[system][snapshot] -> float
    :param sample_size: int
    :param final_fes: dict, {ref_ligand (str): {ligand (str): fe (float)}}
    :return (pearson_r, rmse): (float, float)
    """
    assert ref_ligand in ref_ligands, ref_ligand + " not in " + ref_ligands

    rand_snapshots = np.random.choice(weights[ref_ligand].keys(), size=sample_size, replace=True)

    fes = {}
    for ligand in target_ligands:
        group = ligand[:-3]
        code = ligand[-3:]
        fe_cal = RelBFEWithoutCV(algdock_score_dir, group, code, weights, ref_ligands, yank_interaction_energies)
        fes[ligand] = fe_cal.cal_exp_mean_for_one_ref_ligand(FF, rand_snapshots, ref_ligand)

    pearson_r, rmse = _pearson_r_rmse(final_fes[ref_ligand], fes)
    return pearson_r, rmse


def _r_rmse_one_ref_ligand_a_random_sample_of_snapshot_with_cv_3a(algdock_score_dir, target_ligands,
                                                      ref_ligand, ref_ligands,
                                                      FF, weights, yank_interaction_energies,
                                                      sample_size, final_fes):
    """
    :param algdock_score_dir: str
    :param target_ligands: list of str
    :param ref_ligand: str
    :param FF: str
    :param ref_ligands: list of str
    :param weights: dict,
                        weights[ref_ligand_name][snapshot] -> float
                        weights["systems"][ref_ligand_name] -> float
    :param yank_interaction_energies: dict, yank_interaction_energies[system][snapshot] -> float
    :param sample_size: int
    :param final_fes: dict, {ref_ligand (str): {ligand (str): fe (float)}}
    :return (pearson_r, rmse): (float, float)
    """
    assert ref_ligand in ref_ligands, ref_ligand + " not in " + ref_ligands

    rand_snapshots = np.random.choice(weights[ref_ligand].keys(), size=sample_size, replace=True)

    fes = {}
    for ligand in target_ligands:
        _, _, _, _, fe = relative_bfe_with_cv_using_exp_mean_method_3a(rand_snapshots, algdock_score_dir,
                                                                       ligand, ref_ligand,
                                                                       weights, yank_interaction_energies, FF,
                                                                       remove_outliers_g_h=False,
                                                                       subtract_self=False,
                                                                       flip_sign_c=True,
                                                                       verbose=False)
        fes[ligand] = fe

    pearson_r, rmse = _pearson_r_rmse(final_fes[ref_ligand], fes)
    return pearson_r, rmse


def _bootstrap_r_rmse_one_ref_ligand(algdock_score_dir, target_ligands,
                                     ref_ligand, ref_ligands,
                                     FF, weights, yank_interaction_energies,
                                     sample_size, final_fes, repeats,
                                     method):
    """
    :param algdock_score_dir: str
    :param target_ligands: list of str
    :param ref_ligand: str
    :param ref_ligands: list of str
    :param FF: str
    :param weights:  dict,
                    weights[ref_ligand_name][snapshot] -> float
                    weights["systems"][ref_ligand_name] -> float
    :param yank_interaction_energies: dict, yank_interaction_energies[system][snapshot] -> float
    :param sample_size: int
    :param final_fes: dict, {ref_ligand (str): {ligand (str): fe (float)}}
    :param repeats: int
    :param method: str
    :return (r_mean, r_std, rmse_mean, rmse_std): (float, float, float, float)
    """
    assert method in ["without_cv", "with_cv_3a"], "unrecognized method: " + method
    if method == "without_cv":
        _r_rmse_one_ref_ligand_a_random_sample_of_snapshot = _r_rmse_one_ref_ligand_a_random_sample_of_snapshot_without_cv
    else:
        _r_rmse_one_ref_ligand_a_random_sample_of_snapshot = _r_rmse_one_ref_ligand_a_random_sample_of_snapshot_with_cv_3a

    rs = []
    rmses = []
    for _ in range(repeats):
        r, rmse = _r_rmse_one_ref_ligand_a_random_sample_of_snapshot(algdock_score_dir, target_ligands,
                                                                     ref_ligand, ref_ligands,
                                                                     FF, weights, yank_interaction_energies,
                                                                     sample_size, final_fes)

        if str(r).lower() not in ["-inf", "inf", "nan"]:
            if str(rmse).lower() not in ["-inf", "inf", "nan"]:
                rs.append(r)
                rmses.append(rmse)

    r_mean = np.mean(rs)
    r_std = np.std(rs)
    rmse_mean = np.mean(rmses)
    rmse_std = np.std(rmses)

    return r_mean, r_std, rmse_mean, rmse_std


# sample_sizes = np.array(range(10, 51, 5) + range(60, 91, 10) + [96], dtype=int)
sample_sizes = [int(num) for num in args.sample_sizes.split()]
print("sample_sizes", sample_sizes)

_, _, single_snap_weights, _, _ = load_mbar_weights()

ref_ligands = [ligand for ligand in single_snap_weights.keys() if ligand != "systems"]
print("ref_ligands:", ref_ligands)
target_ligands = YANK_LIGANDS.keys()
print("target_ligands:", target_ligands)

yank_interaction_energies = load_interaction_energies(path=args.interaction_energies_dir)

final_fes = _load_final_fe(args.final_results_dir, ref_ligands, args.weight_scheme, args.combining_rule, args.FF)

for ref_ligand in ref_ligands:
    print("Processing ref ligand:", ref_ligand)
    out_dir = ref_ligand
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    out_file = open(os.path.join(out_dir, "r_rmse.dat"), "w")
    out_file.write("# sample_size     r          r_std          rmse          rmse_std\n")

    for sample_size in sample_sizes:
        print(sample_size)
        r, r_std, rmse, rmse_std = _bootstrap_r_rmse_one_ref_ligand(args.algdock_score_dir, target_ligands,
                                                                    ref_ligand, ref_ligands,
                                                                    args.FF, single_snap_weights, yank_interaction_energies,
                                                                    sample_size, final_fes, args.bootstrap_repeats,
                                                                    args.method)

        out_file.write("%10d %15.10f %15.10f %15.10f %15.10f\n" % (sample_size, r, r_std, rmse, rmse_std))
    out_file.close()

print("DONE")
