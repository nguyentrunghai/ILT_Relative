"""
run convergence with respect to final for relative binding free energies without control variates
"""

import argparse

from _process_yank_outputs import load_interaction_energies

parser = argparse.ArgumentParser()
parser.add_argument("--algdock_score_dir", type=str,
                    default="Bing_Calculations/Correct_Pro_BornRadii/Concatened_Scores/version2016_May_AlGDock/OBC2")
parser.add_argument("--interaction_energies_dir",   type=str,
                    default="Relative_Binding_FE/OpenMM_OBC2_interaction_energies_for_576_algdock_snapshots")

parser.add_argument("--FF", type=str, default="OpenMM_OBC2_MBAR")
parser.add_argument("--weight_scheme", type=str, default="__equal_sys__single_weight")

parser.add_argument("--final_results_dir", type=str,
                    default="Relative_Binding_FE/Relative_FE_Est_1/all96")

args = parser.parse_args()

yank_interaction_energies = load_interaction_energies(path=args.interaction_energies_dir)
