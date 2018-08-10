
import argparse
import glob
import os
import sys
from collections import defaultdict

import mdtraj as md

sys.path.append("/home/tnguye46/T4_Lysozyme/scripts")
from _pca import cal_pca_mult_trajs 

parser = argparse.ArgumentParser()

parser.add_argument( "--holo_traj_glob_matching",  type=str, default="/home/tnguye46/T4_Lysozyme/Yank_snapshots/all_snapshots/heavy_atoms_near_bsite_stride1/*/bound_state0_*.pdb")

parser.add_argument( "--apo_traj_glob_matching",  type=str, default="/home/tnguye46/T4_Lysozyme/Yank_snapshots/all_snapshots/heavy_atoms_near_bsite_stride5/*/unbound_state15_*.pdb")

parser.add_argument( "--apo_out",  type=str, default="apo.dat")
parser.add_argument( "--holo_out_prefix",  type=str, default="holo")

args = parser.parse_args()

def _load_pdbs_into_one_mdtrj(pdb_files):
    assert isinstance(pdb_files, list), "pdb_files must be a list"
    traj = md.load_pdb(pdb_files[0])

    traj = traj.join( [ md.load_pdb(pdb_file) for pdb_file in pdb_files[1:] ], check_topology=True)
    return traj

def _write_pcs(pcs, name, out):
    with open(out, "w") as handle:
        handle.write("# " + name + "\n")
        handle.write("# pc1             pc2\n")
        for pc1, pc2 in pcs:
            handle.write("%20.10f  %20.10f\n"%(pc1, pc2) )
    return None

holo_traj_files = glob.glob(args.holo_traj_glob_matching)
apo_traj_files = glob.glob(args.apo_traj_glob_matching)

holo_traj_files_with_ligands_as_keys = defaultdict(list)
for file_name in holo_traj_files:
    ligand = os.path.dirname(file_name)
    ligand = os.path.basename(ligand)
    holo_traj_files_with_ligands_as_keys[ligand].append(file_name)

combined_apo_md_traj = _load_pdbs_into_one_mdtrj(apo_traj_files)
ligands = holo_traj_files_with_ligands_as_keys.keys()
holo_md_trajs = [_load_pdbs_into_one_mdtrj( holo_traj_files_with_ligands_as_keys[ligand] ) for ligand in ligands]

all_pcs = cal_pca_mult_trajs( [combined_apo_md_traj] + holo_md_trajs, n_components=2)

_write_pcs(all_pcs[0], "apo", args.apo_out)

for ligand, pcs in zip(ligands, all_pcs[1:]):
    out = args.holo_out_prefix + "_" + ligand + ".dat"
    _write_pcs(pcs, ligand, out)

print "DONE"
