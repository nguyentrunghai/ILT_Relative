
import sys
import os
import argparse

import numpy as np

from _averaging_absolute import MultiStruScores 

from _algdock import SIX_YANK_SYSTEMS, FOUR_ACTIVE_YANK_SYSTEMS


from _process_yank_outputs import load_interaction_energies 
from load_mbar_weights_apo_OBC2 import load_mbar_weights

sys.path.append("/home/tnguye46/FFT_T4/scripts")
from _yank_algdock_fft_scores import load_scores

TEMPERATURE = 300.0
KB = 8.3144621E-3/4.184  # kcal/mol/K
BETA = 1. / TEMPERATURE / KB
V_0 = 1661.

def site_correction( r0 ):
    V_binding = 4. / 3. * np.pi * ( r0 ** 3 )
    return -TEMPERATURE * KB * np.log( V_binding / V_0 / 8 / np.pi**2 )

SITE_CORRECTION = site_correction( 5.0 )

def pearson_r_and_rmse(yank_fe, algdock_fe):
    ligands = set(yank_fe.keys()).intersection(algdock_fe.keys())
    ligands = [l for l in ligands if str(algdock_fe[l]) != "inf" ]
    
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

def absolute_fes_using_random_sets_of_receptor_snapshots(score_dir, ligands, weights, six_yank_systems, 
                                                            FF, sample_size, repeats,
                                                            get_min_across_systems=False):
    """
    return a list of dic relative_fes[ligand] -> float
    """
    all_snapshots = get_all_snapshots(six_yank_systems, weights)
    list_of_ramdom_snapshots = [np.random.choice(all_snapshots, size=sample_size, replace=True) for _ in range(repeats)]
    list_of_ramdom_snapshots = [list(ramdom_snapshots) for ramdom_snapshots in list_of_ramdom_snapshots]
    
    absolute_fes = {}
    absolute_fes = [{} for _ in range(repeats)]
    for ligand in ligands:
        #print ligand
        ligand_group   = ligand[:-3]
        ligand_3l_code = ligand[-3:] 
        fe_cal = MultiStruScores(score_dir, ligand_group, ligand_3l_code, weights, six_yank_systems)

        for i, ramdom_snapshots in enumerate(list_of_ramdom_snapshots):

            if get_min_across_systems:
                absolute_fes[i][ligand] = fe_cal.get_min_exp_mean_for_snapshots(FF, ramdom_snapshots) + SITE_CORRECTION
            else:
                absolute_fes[i][ligand] = fe_cal.get_exp_mean_for_snapshots(FF, ramdom_snapshots) + SITE_CORRECTION

        del(fe_cal)

    return absolute_fes

# ----------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--scores_dir",                 type=str, default = "/home/tnguye46/T4_Lysozyme/Bing_Calculations/Correct_Pro_BornRadii/Concatened_Scores/version2016_May_AlGDock/OBC2")

parser.add_argument( "--FF",                    type=str, default = "OpenMM_OBC2_MBAR")
parser.add_argument( "--exclude_ligands",       type=list, default = [])
parser.add_argument( "--yank_results",          type=str, default = "/home/tnguye46/T4_Lysozyme/Yank/yank_results.dat")

parser.add_argument( "--bootstrap_repeats",     type=int, default = 500)

parser.add_argument( "--get_min_across_systems",     type=bool, default = False)
parser.add_argument( "--yank_systems",               type=str, default = "read")

parser.add_argument( "--out_prefix",          type=str, default = "convergence_absolute_average_across_systems")

args = parser.parse_args()

#sample_sizes = np.array( [5] + range(10, 101, 10) + range(150, 551, 50) + [576], dtype=int)

sample_sizes = np.array( range(5, 51, 5) + range(60, 91, 10) + [96], dtype=int)
print sample_sizes

block_weights, state_weights, single_snap_weights, stru_group_weights_equal_sys, stru_group_weights_ub_weighted = load_mbar_weights()

for system in state_weights["systems"]:
    state_weights["systems"][system] = 1.

print state_weights["systems"]

yank_fes, yank_fe_stds = load_scores(args.yank_results, 0, 1, 2, args.exclude_ligands)

ligands = yank_fes.keys()

pearson_r_out_handle = open(args.out_prefix + "_R.dat", "w")
rmse_out_handle      = open(args.out_prefix + "_RMSE.dat", "w")

if args.yank_systems == "all":
    yank_systems = SIX_YANK_SYSTEMS
elif args.yank_systems == "actives":
    yank_systems = FOUR_ACTIVE_YANK_SYSTEMS 
elif args.yank_systems == "read":
    yank_systems = [ os.path.basename( os.getcwd() ) ]
else:
    raise Exception("Unknown yank systems")

print yank_systems

for sample_size in sample_sizes:

    absolute_fes = absolute_fes_using_random_sets_of_receptor_snapshots(args.scores_dir, ligands, state_weights, 
                                                                        yank_systems, args.FF, sample_size, args.bootstrap_repeats,
                                                                        get_min_across_systems=args.get_min_across_systems)
    rs = []
    rmses = []
    for absolute_fe in absolute_fes:
        r, rmse = pearson_r_and_rmse(yank_fes, absolute_fe)
        rs.append(r)
        rmses.append(rmse)

    rs = [r for r in rs if str(r).lower() not in ["inf", "nan"] ]
    rmses = [rmse for rmse in rmses if str(rmse).lower() not in ["inf", "nan"] ]
    print sample_size, np.mean(rs), np.mean(rmses)

    pearson_r_out_handle.write("%10d %20.10f %20.10f\n" %(sample_size, np.mean(rs), np.std(rs)) )
    rmse_out_handle.write("%10d %20.10f %20.10f\n" %(sample_size, np.mean(rmses), np.std(rmses)) )

print "Done"


