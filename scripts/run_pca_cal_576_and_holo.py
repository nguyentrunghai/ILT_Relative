
import argparse
import glob
import os
import sys

import mdtraj as md

from _algdock import SIX_YANK_SYSTEMS, load_algdock_snapshots_for_each_of_six_yank_systems
from load_mbar_weights_holo_OBC2 import load_mbar_weights

sys.path.append("/home/tnguye46/T4_Lysozyme/scripts")
from _pca import cal_pca_use_ref

parser = argparse.ArgumentParser()
parser.add_argument( "--algdock_576_snap_dir",       type=str, default = "/home/tnguye46/T4_Lysozyme/Yank_snapshots/576_rec_snapshots/heavy_atoms_near_bsite")
parser.add_argument( "--all_24_yank_bound_snap",       type=str, default = "/home/tnguye46/T4_Lysozyme/Yank_snapshots/all_snapshots/heavy_atoms_near_bsite/*/bound_state0*pdb")

parser.add_argument( "--ref_out",       type=str, default ="all_24_holo.dat")

args = parser.parse_args()

algdock_snapshots_for_each_of_six_yank_systems = load_algdock_snapshots_for_each_of_six_yank_systems()
# take single_snap_weights
block_weights, state_weights, single_snap_weights, stru_group_weights_equal_sys, stru_group_weights_ub_weighted = load_mbar_weights() 

def _load_pdbs_into_one_mdtrj(pdb_files):
    assert isinstance(pdb_files, list), "pdb_files must be a list"
    traj = md.load_pdb(pdb_files[0])

    traj = traj.join( [ md.load_pdb(pdb_file) for pdb_file in pdb_files[1:] ], check_topology=True)
    return traj

ref_pdb_files = glob.glob(args.all_24_yank_bound_snap)
ref_traj = _load_pdbs_into_one_mdtrj(ref_pdb_files)

trajs = []
for yank_system in SIX_YANK_SYSTEMS:
    print yank_system
    snapshots = algdock_snapshots_for_each_of_six_yank_systems[yank_system]

    pdb_files = [ os.path.join(args.algdock_576_snap_dir, snapshot + ".pdb") for snapshot in snapshots ]

    traj = _load_pdbs_into_one_mdtrj(pdb_files)

    trajs.append(traj)

bound_pcs, algdock_pcs = cal_pca_use_ref([ref_traj], trajs, 2)

with open(args.ref_out, "w") as handle:
    handle.write("# all 24 bound\n")
    handle.write("# pc1         pc2\n")
    for pcs in bound_pcs:
        for pc1, pc2 in pcs: 
            handle.write("%20.10f  %20.10f\n"%(pc1, pc2) )

for yank_system, pcs in zip(SIX_YANK_SYSTEMS, algdock_pcs):
    
    out = yank_system + "_pc.dat"

    with open(out, "w") as handle:
        handle.write("# " + yank_system + "\n")
        handle.write("# pc1          pc2        holo weight\n")

        for i, snapshot in enumerate( algdock_snapshots_for_each_of_six_yank_systems[yank_system] ):
            pc1, pc2 = pcs[i]
            weight = single_snap_weights[yank_system][snapshot]
            handle.write("%20.10f  %20.10f %20.10e\n"%(pc1, pc2, weight) )


