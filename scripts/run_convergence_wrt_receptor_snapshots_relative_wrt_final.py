
import os
import sys
import argparse

import numpy as np

from _relative_estimators import RelBFEWithoutCV

from _algdock import SIX_YANK_SYSTEMS
from _process_yank_outputs import load_interaction_energies 
from load_mbar_weights_holo_OBC2 import load_mbar_weights

sys.path.append("/home/tnguye46/FFT_T4/scripts")
from _yank import load_scores


def pearson_r_and_rmse(yank_fe, algdock_fe):
    """
    :param yank_fe: dict, ligand (str) -> free energy (float)
    :param algdock_fe: dict, ligand (str) --> free energy (float)
    :return (r, rmse): (float, float)
    """
    ligands = set(yank_fe.keys()).intersection(algdock_fe.keys())
    ligands = [ligand for ligand in ligands if str(algdock_fe[ligand]).lower() not in ["inf", "nan"]]
    
    xd = np.array([yank_fe[ligand] for ligand in ligands], dtype=float)
    yd = np.array([algdock_fe[ligand] for ligand in ligands], dtype=float)

    r = np.corrcoef( [xd, yd] )[0][-1]
    rmse = ((xd - yd)**2).mean()
    rmse = np.sqrt(rmse)

    return r, rmse


def get_all_snapshots(six_yank_systems, weights):
    """
    :param six_yank_systems: list of str
    :param weights: dict, ligand (str) --> {snapshot (str) : weight (float)}
    :return:
    """
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

def relative_fes_final_results(score_dir, ligands, weights, six_yank_systems, yank_interaction_energies, FF):
    all_snapshots = get_all_snapshots(six_yank_systems, weights)

    relative_fes = {}
    for ligand in ligands:
        ligand_group   = ligand[:-3]
        ligand_3l_code = ligand[-3:] 
        fe_cal = RelBFEWithoutCV(score_dir, ligand_group, ligand_3l_code, weights, six_yank_systems, yank_interaction_energies)

        relative_fes[ligand] = fe_cal.cal_exp_mean_separate_for_each_system(FF, all_snapshots)
    return relative_fes

def relative_to_absolute(relative_fes, yank_fes):
    """
    this take care of non zero self relative
    """
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

parser.add_argument( "--final_results_dir",     type=str, default = "/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/Relative_FE_Est_1/OBC2")
parser.add_argument( "--weight_scheme",         type=str, default = "__equal_sys__single_weight")
parser.add_argument( "--file_to_take_all_ligands",         type=str, default = "/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/Relative_FE_Est_1/OBC2/1-methylpyrrole.A__AAA__equal_sys__single_weight/ExpMean/OpenMM_OBC2_MBAR.score")

parser.add_argument( "--bootstrap_repeats",     type=int, default = 500)

parser.add_argument( "--yank_systems",               type=str, default = "all") 

parser.add_argument( "--out_prefix",          type=str, default = "convergence_relative_est_1")
#parser.add_argument( "--out_prefix",          type=str, default = "convergence_relative_not_convert_2_abs_est_1")

args = parser.parse_args()

sample_sizes = np.array( [5] + range(10, 101, 10) + range(150, 551, 50) + [576], dtype=int)
#sample_sizes = np.array( range(5, 51, 5) + range(60, 91, 10) + [96], dtype=int)
print sample_sizes

block_weights, state_weights, single_snap_weights, stru_group_weights_equal_sys, stru_group_weights_ub_weighted = load_mbar_weights()
yank_interaction_energies = load_interaction_energies(path=args.interaction_energies_dir) 

yank_fes, yank_fe_stds = load_scores(args.yank_results, 0, 1, 2, args.exclude_ligands)

if args.yank_systems == "all":
    yank_systems = SIX_YANK_SYSTEMS
elif args.yank_systems == "actives":
    yank_systems = FOUR_ACTIVE_YANK_SYSTEMS 
elif args.yank_systems == "read":
    yank_systems = [ os.path.basename( os.getcwd() ) ]
else:
    raise Exception("Unknown yank systems")

if len(yank_systems) > 1:
    if args.yank_systems != "all":
        raise Exception("More than one yank systems but not all 6")
    else:
        print "Loading " + args.file_to_take_all_ligands
        dummy_results, _ = load_scores(args.file_to_take_all_ligands, 0, 1, 2, args.exclude_ligands)
        ligands = dummy_results.keys()

        final_results = relative_fes_final_results(args.scores_dir, ligands, single_snap_weights, yank_systems, yank_interaction_energies, args.FF) 

elif len(yank_systems) == 1:
    ref_ligand = yank_systems[0]
    final_results_file = os.path.join(args.final_results_dir, ref_ligand + args.weight_scheme, "ExpMean", args.FF + ".score")
    print "Loading " + final_results_file

    un_normalized_rel_fes, _ = load_scores(final_results_file, 0, 1, 2, args.exclude_ligands)
    final_results = { ligand : {ref_ligand : un_normalized_rel_fes[ligand]} for ligand in un_normalized_rel_fes}
    

else:
    raise Exception("something is wrong")

# TODO
final_results = relative_to_absolute( final_results,  yank_fes)
final_results = min_absolute_fes(final_results)
ligands = final_results.keys()

pearson_r_out_handle = open(args.out_prefix + "_R.dat", "w")
rmse_out_handle      = open(args.out_prefix + "_RMSE.dat", "w")

for sample_size in sample_sizes:
    relative_fes = relative_fes_using_random_sets_of_receptor_snapshots(args.scores_dir, ligands, single_snap_weights, yank_systems, 
                                                                            yank_interaction_energies, args.FF, sample_size, args.bootstrap_repeats)
    rs = []
    rmses = []
    for relative_fe in relative_fes:
        # TODO
        abs_fe     = relative_to_absolute(relative_fe, yank_fes)
        min_abs_fe = min_absolute_fes(abs_fe)
        #min_abs_fe = min_absolute_fes(relative_fe)

        r, rmse = pearson_r_and_rmse(final_results, min_abs_fe)
        rs.append(r)
        rmses.append(rmse)


    rs = [r for r in rs if str(r).lower() not in ["inf", "nan"] ]
    rmses = [rmse for rmse in rmses if str(rmse).lower() not in ["inf", "nan"] ]
    print sample_size, np.mean(rs), np.mean(rmses)

    pearson_r_out_handle.write("%10d %20.10f %20.10f\n" %(sample_size, np.mean(rs), np.std(rs)) )
    rmse_out_handle.write("%10d %20.10f %20.10f\n" %(sample_size, np.mean(rmses), np.std(rmses)) )

print "Done"

