#
import os
import sys
import argparse
import glob
import numpy as np
import copy
np.seterr(over='raise') # raise exception if overload

from _algdock import load_algdock_snapshots_for_each_of_six_yank_systems

TEMPERATURE = 300.                                                                                                                                                       
KB = 8.3144621E-3/4.184  # kcal/mol/K
BETA = 1. / TEMPERATURE / KB

dock6_dir   = "dock6"
algdock_dir = "AlGDock/dock"
#algdock_dir = "AlGDock"

EXCLUDE_FFS = ["receptor_OpenMM_Gas", "receptor_OpenMM_OBC2", "Theta_1OpenMM_Gas", "Theta_1OpenMM_OBC2", "Theta_RLOpenMM_Gas", "Theta_RLOpenMM_OBC2", 
        "Theta_1sander_Gas", "Theta_1sander_PBSA", "Theta_RLsander_Gas", "Theta_RLsander_PBSA", "receptor_sander_Gas", "receptor_sander_PBSA", 
        "grid_MBAR", "MBAR", "OpenMM_Gas_MBAR_c2", "OpenMM_Gas_inverse_FEP", "OpenMM_OBC2_inverse_FEP", "OpenMM_OBC2_MBAR_c2"]

class MultiStruScores:
    """
    load and averaging scores for a ligand over multiple receptor's snapshots
    """
    def __init__( self, score_dir, ligand_group, ligand_3l_code, weights, yank_systems, exclude_ffs=EXCLUDE_FFS, 
                repeats=100):
        """
        weights[system][snapshot]
        yank_systems:   list of str
        """
        self._identification = ligand_group + ligand_3l_code
        self._dock_dir = os.path.join( score_dir, dock6_dir, ligand_group, ligand_3l_code )
        self._algdock_dir = os.path.join( score_dir, algdock_dir, ligand_group, ligand_3l_code )
        
        self._FFs = ["dock6"]
        FFs = glob.glob( os.path.join(self._algdock_dir, "*.score") )
        FFs = [ os.path.basename(FF)[:-6] for FF in FFs ]
        FFs = [ff for ff in FFs if ff not in exclude_ffs]
        self._FFs.extend( FFs )
        
        self._scores = {}
        self._load_dock6()
        self._load_algdock()
         
        self._weights = weights
        self._yank_systems = yank_systems
        self._repeats = repeats
        self._considered_snapshots = self._get_snapshots()
        self._allowed_snashots = self._get_snapshots_in_scores_and_systems()

    def _load_dock6(self):
        in_file = open( os.path.join(self._dock_dir, "dock6.score"), "r")
        entries = {}
        for line in in_file:
            words = line.split()
            snap_id, value = words[0], words[1]
            
            if value.lower() != "nan":
                entries[snap_id] = np.float(value) / TEMPERATURE / KB

        self._scores["dock6"] = entries
        in_file.close()
        return None
    
    def _load_algdock(self):
        for FF in self._FFs:
            if FF != "dock6":
                in_file = open( os.path.join(self._algdock_dir, FF + ".score"), "r" )
                entries = {}
                for line in in_file:
                    words = line.split()
                    snap_id, value = words[0], words[1]

                    if value.lower() != "nan":
                        entries[snap_id] = np.float(value)

                self._scores[FF] = entries
                in_file.close()
        return None

    def _get_snapshots(self):
        snapshots = []
        for system in self._yank_systems:
            for snapshot in self._weights[system].keys():
                snapshots.append(snapshot)
        return snapshots

    def _get_snapshots_in_scores_and_systems(self):
        snapshots = {}
        for FF in self._FFs:
            snapshots[FF] = {}
            snapshots_have_scores = self._scores[FF].keys()

            for system in self._yank_systems:
                snapshots_in_system = self._weights[system].keys()
                snapshots[FF][system] = set(snapshots_have_scores).intersection(snapshots_in_system)
        return snapshots
    
    def get_id(self):
        return self._identification
    
    def get_FFs(self):
        return self._FFs
    
    def _cal_exp_mean(self, snapshots):
        averages = {}
        for FF in self._FFs:
            sys_mean = 0.
            for system in self._yank_systems:
                a = 0.
                w = 0.
                for snapshot in snapshots:
                    if snapshot in self._allowed_snashots[FF][system]:
                        try:
                            a += np.exp( -1.0 * self._scores[FF][snapshot] ) * self._weights[system][snapshot]
                        except FloatingPointError:
                            print "overflow for " + self._identification + " " + snapshot + " " + FF
                        else:
                            w += self._weights[system][snapshot]
                if w != 0:
                    a = a / w
                sys_mean += a * self._weights["systems"][system]

            sys_mean = sys_mean / sum( [self._weights["systems"][system] for system in self._yank_systems] )
            averages[FF] = (-1./BETA) * np.log(sys_mean)
        return averages

    def get_exp_mean(self):
        return self._cal_exp_mean(self._considered_snapshots)

    def get_exp_mean_std(self):
        fes = {FF:[] for FF in self._FFs}
        for repeat in range(self._repeats):
            snapshots = np.random.choice(self._considered_snapshots, size=len(self._considered_snapshots), replace=True)
            fe = self._cal_exp_mean(snapshots)

            for FF in self._FFs:
                if fe[FF] not in [np.inf, -np.inf, np.nan]:
                    fes[FF].append(fe[FF])
        std = {}
        for FF in self._FFs:
            if len(fes[FF]) > 0:
                std[FF] = np.std(fes[FF])
            else:
                std[FF] = 0.
        return std
    
    def _cal_mean(self, snapshots):
        averages = {}
        for FF in self._FFs:
            sys_mean = 0.
            for system in self._yank_systems:
                a = 0.
                w = 0.
                for snapshot in snapshots:
                    if snapshot in self._allowed_snashots[FF][system]:
                        if self._scores[FF][snapshot] != np.inf:
                            a += self._scores[FF][snapshot] * self._weights[system][snapshot]
                            w += self._weights[system][snapshot]
                if w != 0:
                    a = a / w
                sys_mean += a * self._weights["systems"][system]

            sys_mean = sys_mean / np.sum( [self._weights["systems"][system] for system in self._yank_systems] )
            averages[FF] = sys_mean * TEMPERATURE * KB
        return averages

    def get_mean(self):
        return self._cal_mean(self._considered_snapshots)

    def get_mean_std(self):
        fes = {FF:[] for FF in self._FFs}
        for repeat in range(self._repeats):
            snapshots = np.random.choice(self._considered_snapshots, size=len(self._considered_snapshots), replace=True)
            fe = self._cal_mean(snapshots)

            for FF in self._FFs:
                if fe[FF] not in [np.inf, -np.inf, np.nan]:
                    fes[FF].append(fe[FF])
        std = {}
        for FF in self._FFs:
            if len(fes[FF]) > 0:
                std[FF] = np.std(fes[FF])
            else:
                std[FF] = 0.
        return std 
    
    def _cal_min(self, snapshots):
        averages = {}
        for FF in self._FFs:
            a = []
            for system in self._yank_systems:
                a += [ self._scores[FF][snapshot] for snapshot in snapshots if snapshot in self._allowed_snashots[FF][system] ]

            if len(a) > 0:
                averages[FF] = np.array(a).min() * TEMPERATURE * KB
            else:
                averages[FF] = np.inf
        return averages

    def get_min(self):
        return self._cal_min(self._considered_snapshots)

    def get_min_std(self):
        fes = {FF:[] for FF in self._FFs}
        for repeat in range(self._repeats):
            snapshots = np.random.choice(self._considered_snapshots, size=len(self._considered_snapshots), replace=True)
            fe = self._cal_min(snapshots)

            for FF in self._FFs:
                if fe[FF] not in [np.inf, -np.inf, np.nan]:
                    fes[FF].append(fe[FF])
        std = {}
        for FF in self._FFs:
            if len(fes[FF]) > 0:
                std[FF] = np.std(fes[FF])
            else:
                std[FF] = 0.
        return std

    def check_extreme_low(self, cutoff = -100.):
        for FF in self._FFs:
            for system in self._yank_systems:
                for snapshot in self._weights[system].keys():
                    if snapshot in self._allowed_snashots[FF][system]:
                        if self._scores[FF][snapshot] != np.inf:
                            if self._scores[FF][snapshot] < cutoff:
                                print "Extreme low: " + self._identification + " " + snapshot + " " + FF + "  %20.10f"%self._scores[FF][snapshot]
        return None
#----------------------------
#------------------------------------------------------

def averaging( yank_systems, result_dir, weights ):
    #
    if not os.path.isdir(result_dir):
        os.system("mkdir " + result_dir)

    out_files = {}
    for rule in combining_rules:
        if not os.path.isdir( os.path.join(result_dir, rule) ):
            os.system("mkdir " + os.path.join(result_dir, rule) )
        out_files[rule] = {}
        for FF in FFs:
            out_files[rule][FF] = open( os.path.join(result_dir, rule, FF + '.score'), 'w' )

    for group in ligand_3l_codes.keys():
        for code in ligand_3l_codes[group]:
            scores = MultiStruScores(args.scores_dir, group, code,  weights, yank_systems )
            scores.check_extreme_low()
            for rule in combining_rules:
                if rule == 'Mean':
                    averages     = scores.get_mean()
                    standard_dev = scores.get_mean_std()

                elif rule == 'Min':
                    averages     = scores.get_min()
                    standard_dev = scores.get_min_std()

                elif rule == 'ExpMean':
                    averages     = scores.get_exp_mean()
                    standard_dev = scores.get_exp_mean_std()
                else:
                    raise RuntimeError('unknown combining rule')

                id = scores.get_id()
                for FF in averages.keys():
                    out_files[rule][FF].write("%s   %20.10f %20.10f\n" %(id, averages[FF], standard_dev[FF]) )

    for rule in combining_rules:
        for FF in FFs:
            out_files[rule][FF].close()

    return None

def equalize_system_weights( original_weights ):
    new_weights = copy.deepcopy( original_weights )
    for system in new_weights["systems"].keys():
        new_weights["systems"][system] = 1.0
    return new_weights


def take_6_apo(original_weights):
    new_weights = equalize_system_weights(original_weights)
    ordered_snapshots = load_algdock_snapshots_for_each_of_six_yank_systems() 

    for system in new_weights["systems"].keys():
        for i, snapshot in enumerate( ordered_snapshots[system] ):
            if i < len(ordered_snapshots[system]) - 6:
                new_weights[system][snapshot] = 0.
            else:
                new_weights[system][snapshot] = 1.
            
    return new_weights


def take_12_near_apo(original_weights):
    new_weights = equalize_system_weights(original_weights)
    ordered_snapshots = load_algdock_snapshots_for_each_of_six_yank_systems()

    for system in new_weights["systems"].keys():
        for i, snapshot in enumerate( ordered_snapshots[system] ):
            if i < len(ordered_snapshots[system]) - 12:
                new_weights[system][snapshot] = 0.
    return new_weights


def take_24_near_apo(original_weights):
    new_weights = equalize_system_weights(original_weights)
    ordered_snapshots = load_algdock_snapshots_for_each_of_six_yank_systems()

    for system in new_weights["systems"].keys():
        for i, snapshot in enumerate( ordered_snapshots[system] ):
            if i < len(ordered_snapshots[system]) - 24:
                new_weights[system][snapshot] = 0.
    return new_weights


#---------------------

parser = argparse.ArgumentParser()
parser.add_argument("--scores_dir", type=str, default = "/home/tnguye46/T4_Lysozyme/Bing_Calculations/Correct_Pro_BornRadii/Concatened_Scores/version2016_May_AlGDock/OBC2")

parser.add_argument("--which_snapshots_to_take", type=str, default="all96")    # "all96", "6apo", "12nearapo", "24nearapo" 
args = parser.parse_args()
#

assert args.which_snapshots_to_take in ["all96", "6apo", "12nearapo", "24nearapo"], "unknown which_snapshots_to_take"
print "use " + args.which_snapshots_to_take

combining_rules = [ 'Mean', 'ExpMean', 'Min' ]
#
ligand_groups = glob.glob( os.path.join(args.scores_dir, dock6_dir, "*") )
ligand_groups = [os.path.basename(dir) for dir in ligand_groups]
#
if os.path.basename(args.scores_dir) == "OBC2":
    print "OBC2 weights"
    from load_mbar_weights_apo_OBC2 import load_mbar_weights
elif os.path.basename(args.scores_dir) == "PBSA":
    print "PBSA weights"
    from load_mbar_weights_apo_PBSA import load_mbar_weights
else:
    raise RuntimeError("unknown phase "+os.path.basename(args.scores_dir))
#
ligand_3l_codes = {}
for group in ligand_groups:
    codes = glob.glob( os.path.join( args.scores_dir, dock6_dir, group, "*") )
    ligand_3l_codes[group] = [os.path.basename(c) for c in codes]



block_weights, state_weights, single_snap_weights, stru_group_weights_equal_sys, stru_group_weights_ub_weighted = load_mbar_weights()

yank_systems = [key for key in block_weights.keys() if key not in ["systems"] ]
print "yank systems ", yank_systems


state_weights_equal_systems = equalize_system_weights(state_weights)

print ligand_groups[0], ligand_3l_codes[ligand_groups[0]][0]
FFs = MultiStruScores(args.scores_dir, ligand_groups[0], ligand_3l_codes[ligand_groups[0]][0], block_weights, yank_systems).get_FFs()
print FFs

if args.which_snapshots_to_take == "all96":
    use_state_weights  = equalize_system_weights(state_weights)

elif args.which_snapshots_to_take == "6apo":
    use_state_weights  = take_6_apo(state_weights)

elif args.which_snapshots_to_take == "12nearapo":
    use_state_weights  = take_12_near_apo(state_weights)

elif args.which_snapshots_to_take == "24nearapo":
    use_state_weights  = take_24_near_apo(state_weights)

averaging(yank_systems, "all_snap__equal_sys__state_weight", use_state_weights)


