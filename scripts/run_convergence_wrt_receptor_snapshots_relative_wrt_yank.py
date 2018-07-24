
import os
import sys
import argparse

import numpy as np

from _averaging_relative_estimator import MultiStruScores

from _algdock import SIX_YANK_SYSTEMS
from _process_yank_outputs import load_interaction_energies 
from load_mbar_weights_holo_OBC2 import load_mbar_weights

sys.path.append("/home/tnguye46/FFT_T4/scripts")
from _yank_algdock_fft_scores import load_scores

def pearson_r_and_rmse(yank_fe, algdock_fe):
    ligands = set(yank_fe.keys()).intersection(algdock_fe.keys())
    ligands = [l for l in ligands if str(algdock_fe[l]).lower() not in ["inf", "nan"] ]
    
    xd = np.array([yank_fe[l] for l in ligands], dtype=float)
    yd = np.array([algdock_fe[l] for l in ligands], dtype=float)

    r = np.corrcoef( [xd, yd] )[0][-1]

    rmse = ((xd - yd)**2).mean()
    rmse = np.sqrt(rmse)
    return r, rmse

def get_all_snapshots(six_yank_systems, weights):
    snapshots = []
    for system in six_yank_systems:
        snapshots += weights[system].keys()
    return snapshots

def relative_fes_using_random_sets_of_receptor_snapshots(score_dir, ligands, weights, six_yank_systems, yank_interaction_energies, 
                                                            FF, sample_size, repeats):
    """
    return a list of dic relative_fes[ligand][relative_to_ligand] -> float
    """
    all_snapshots = get_all_snapshots(six_yank_systems, weights)
    list_of_ramdom_snapshots = [np.random.choice(all_snapshots, size=sample_size, replace=True) for _ in range(repeats)]
    list_of_ramdom_snapshots = [list(ramdom_snapshots) for ramdom_snapshots in list_of_ramdom_snapshots]
    
    relative_fes = {}
    relative_fes = [{} for _ in range(repeats)]
    for ligand in ligands:
        #print ligand
        ligand_group   = ligand[:-3]
        ligand_3l_code = ligand[-3:] 
        fe_cal = MultiStruScores(score_dir, ligand_group, ligand_3l_code, weights, six_yank_systems, yank_interaction_energies)

        for i, ramdom_snapshots in enumerate(list_of_ramdom_snapshots):

            #print i
            relative_fes[i][ligand] = fe_cal.cal_exp_mean_separate_for_each_system(FF, ramdom_snapshots)

        del(fe_cal)

    return relative_fes

def relative_to_absolute(relative_fes, yank_fes):
    abs_fes = {}
    for ligand in relative_fes:

        abs_fes[ligand] = {}
        for ref_ligand in relative_fes[ligand]:

            if relative_fes[ref_ligand][ref_ligand] == np.inf:
                val = np.inf
            else:
                val = relative_fes[ligand][ref_ligand] - relative_fes[ref_ligand][ref_ligand] + yank_fes[ref_ligand]

            abs_fes[ligand][ref_ligand] = val
    return abs_fes

def min_absolute_fes(abs_fes):
    min_abs_fes = {}
    for ligand in abs_fes:
        abs_fes[ligand] = np.min( abs_fes[ligand].values() )
    return abs_fes

# ----------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--scores_dir",                 type=str, default = "/home/tnguye46/T4_Lysozyme/Bing_Calculations/Correct_Pro_BornRadii/Concatened_Scores/version2016_May_AlGDock/OBC2")
parser.add_argument("--interaction_energies_dir",   type=str, default = "/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/OpenMM_OBC2_interaction_energies_for_576_algdock_snapshots")

parser.add_argument( "--FF",                    type=str, default = "OpenMM_OBC2_MBAR")
parser.add_argument( "--exclude_ligands",       type=list, default = [])
parser.add_argument( "--yank_results",          type=str, default = "/home/tnguye46/T4_Lysozyme/Yank/yank_results.dat")

parser.add_argument( "--bootstrap_repeats",     type=int, default = 500)

parser.add_argument( "--yank_systems",               type=str, default = "read") 

parser.add_argument( "--out_prefix",          type=str, default = "convergence_relative_est_1")

args = parser.parse_args()

#sample_sizes = np.array( [5] + range(10, 101, 10) + range(150, 551, 50) + [576], dtype=int)
sample_sizes = np.array( range(5, 51, 5) + range(60, 91, 10) + [96], dtype=int)
print sample_sizes

block_weights, state_weights, single_snap_weights, stru_group_weights_equal_sys, stru_group_weights_ub_weighted = load_mbar_weights()
yank_interaction_energies = load_interaction_energies(path=args.interaction_energies_dir) 

yank_fes, yank_fe_stds = load_scores(args.yank_results, 0, 1, 2, args.exclude_ligands)

ligands = yank_fes.keys()

if args.yank_systems == "all":
    yank_systems = SIX_YANK_SYSTEMS
elif args.yank_systems == "actives":
    yank_systems = FOUR_ACTIVE_YANK_SYSTEMS 
elif args.yank_systems == "read":
    yank_systems = [ os.path.basename( os.getcwd() ) ]
else:
    raise Exception("Unknown yank systems")

print yank_systems


pearson_r_out_handle = open(args.out_prefix + "_R.dat", "w")
rmse_out_handle      = open(args.out_prefix + "_RMSE.dat", "w")

for sample_size in sample_sizes:
    relative_fes = relative_fes_using_random_sets_of_receptor_snapshots(args.scores_dir, ligands, single_snap_weights, yank_systems, 
                                                                            yank_interaction_energies, args.FF, sample_size, args.bootstrap_repeats)
    rs = []
    rmses = []
    for relative_fe in relative_fes:
        abs_fe     = relative_to_absolute(relative_fe, yank_fes)
        min_abs_fe = min_absolute_fes(abs_fe)

        r, rmse = pearson_r_and_rmse(yank_fes, min_abs_fe)
        rs.append(r)
        rmses.append(rmse)


    rs = [r for r in rs if str(r).lower() not in ["inf", "nan"] ]
    rmses = [rmse for rmse in rmses if str(rmse).lower() not in ["inf", "nan"] ]
    print sample_size, np.mean(rs), np.mean(rmses)

    pearson_r_out_handle.write("%10d %20.10f %20.10f\n" %(sample_size, np.mean(rs), np.std(rs)) )
    rmse_out_handle.write("%10d %20.10f %20.10f\n" %(sample_size, np.mean(rmses), np.std(rmses)) )

print "Done"

