"""
run convergence with respect to final for relative binding free energies without control variates
"""

import os
import argparse

from _process_yank_outputs import load_interaction_energies
from load_mbar_weights_holo_OBC2 import load_mbar_weights
from _yank import load_scores

parser = argparse.ArgumentParser()
parser.add_argument("--algdock_score_dir", type=str,
                    default="/home/tnguye46/T4_Lysozyme/Bing_Calculations/Correct_Pro_BornRadii/Concatened_Scores/version2016_May_AlGDock/OBC2")
parser.add_argument("--interaction_energies_dir",   type=str,
                    default="/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/OpenMM_OBC2_interaction_energies_for_576_algdock_snapshots")

parser.add_argument("--FF", type=str, default="OpenMM_OBC2_MBAR")
parser.add_argument("--weight_scheme", type=str, default="__equal_sys__single_weight")
parser.add_argument("--combining_rule", type=str, default="ExpMean")

parser.add_argument("--final_results_dir", type=str,
                    default="/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/Relative_FE_Est_1/all96")

args = parser.parse_args()


def _load_final_fe(final_results_dir, ref_ligands, weight_scheme, combining_rule, FF):
    """
    :param final_results_dir: str
    :param ref_ligands: list of str
    :param weight_scheme: str
    :param combining_rule: str
    :param score_file: str
    :return:
    """
    final_fes = {}
    for ref_ligand in ref_ligands:
        fe_file = os.path.join(final_results_dir, ref_ligand+weight_scheme, combining_rule, FF+".score")
        print("Loading file for final fes:", fe_file)
        fes, _ = load_scores(fe_file, 0, 1, 2, exclude_ligands=[])
        fes.pop(ref_ligand)

        final_fes[ref_ligand] = fes
    return final_fes


_, _, single_snap_weights, _, _ = load_mbar_weights()
ref_ligands = [ligand for ligand in single_snap_weights.keys() if ligand != "systems"]
print("ref_ligands", ref_ligands)
yank_interaction_energies = load_interaction_energies(path=args.interaction_energies_dir)

final_fes = _load_final_fe(args.final_results_dir, ref_ligands, args.weight_scheme, args.combining_rule, args.FF)
