
import argparse
import os
import numpy as np

from _process_yank_outputs import read_reduced_energies_for_algdock_snapshots
from _process_yank_outputs import load_receptor_pot_energies
from _process_yank_outputs import augmented_energy_matrix_for_holo_marginal_weights

from _algdock import SIX_YANK_SYSTEMS
from _algdock import load_bpmfs
from _algdock import load_algdock_snapshots_for_each_of_six_yank_systems

from _mbar import mbar_weights_in_state

parser = argparse.ArgumentParser()
parser.add_argument('--yank_dir',       type=str , default="/home/tnguye46/T4_Lysozyme/Yank")
parser.add_argument('--yank_repeats',   type=int, default=3)
parser.add_argument('--yank_nc_file',   type=str, default='complex-implicit.nc')

parser.add_argument('--algdock_dir',   type=str, default='/home/tnguye46/T4_Lysozyme/Bing_Calculations/Correct_Pro_BornRadii/Concatened_Scores/version2016_May_AlGDock/OBC2/AlGDock/dock')
parser.add_argument('--bpmfs_file',   type=str, default='OpenMM_OBC2_MBAR.score')

parser.add_argument('--receptor_pot_energies_dir',  
        type=str , default="/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/OpenMM_OBC2_interaction_energies_for_576_algdock_snapshots")
args = parser.parse_args()

snapshots_for_each_system = load_algdock_snapshots_for_each_of_six_yank_systems()

for yank_ligand in SIX_YANK_SYSTEMS:
    yank_nc_files = [ os.path.join(args.yank_dir, yank_ligand, "Repeat%d"%repeat, "output", args.yank_nc_file)
                            for repeat in range(args.yank_repeats) ]

    u_kln, N_k = read_reduced_energies_for_algdock_snapshots(yank_nc_files)

    receptor_pot_energies_file = os.path.join(args.receptor_pot_energies_dir, yank_ligand+".dat")
    receptor_pot_energies = load_receptor_pot_energies(receptor_pot_energies_file)

    group = yank_ligand[:-3]
    code  = yank_ligand[-3:]
    bpmfs_file = os.path.join(args.algdock_dir, group, code, args.bpmfs_file)
    print "loading " + bpmfs_file
    bpmfs = load_bpmfs(bpmfs_file, exclude_nan=False)

    snapshots = snapshots_for_each_system[yank_ligand]

    bpmfs = [bpmfs[snapshot] for snapshot in snapshots]
    bpmfs = np.array(bpmfs)

    aug_u_kln, aug_N_k = augmented_energy_matrix_for_holo_marginal_weights(u_kln, N_k, bpmfs, receptor_pot_energies)

    weights = mbar_weights_in_state(aug_u_kln, aug_N_k, -1)

    out_file = yank_ligand + ".weig"
    with open(out_file, "w") as handle:
        for snapshot, weight in zip(snapshots, weights):
            handle.write("%10s %20.10e\n" % (snapshot, weight))


