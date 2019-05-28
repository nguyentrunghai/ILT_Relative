"""
"""

import argparse

from _yank import YANK_LIGANDS as target_ligands
from load_mbar_weights_holo_OBC2 import load_mbar_weights
from _process_yank_outputs import load_interaction_energies

parser = argparse.ArgumentParser()

parser.add_argument("--scores_dir", type=str,
                    default="Bing_Calculations/Correct_Pro_BornRadii/Concatened_Scores/version2016_May_AlGDock/OBC2")

parser.add_argument("--interaction_energies_dir", type=str,
                    default="Relative_Binding_FE/OpenMM_OBC2_interaction_energies_for_576_algdock_snapshots")

parser.add_argument("--FF", type=str, default="OpenMM_OBC2_MBAR")

parser.add_argument("--bootstrap_repeats", type=int, default=1000)

args = parser.parse_args()

_, _, single_snap_weights, _, _ = load_mbar_weights()
ref_ligands = [ligand for ligand in single_snap_weights.keys() if ligand != "systems"]
print("ref_ligands", ref_ligands)

yank_interaction_energies = load_interaction_energies(path=args.interaction_energies_dir)

for ref_ligand in ref_ligands:
    snapshots = single_snap_weights[ref_ligand].keys()qui