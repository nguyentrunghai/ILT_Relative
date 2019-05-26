
from __future__ import print_function

import os
import glob
import numpy as np
np.seterr(over='raise') # raise exception if overload

TEMPERATURE = 300.                                                                                                                                                       
KB = 8.3144621E-3/4.184  # kcal/mol/K
BETA = 1. / TEMPERATURE / KB

DOCK6_SUB_DIR = "dock6"
ALGDOCK_SUB_DIR = "AlGDock/dock"

EXCLUDE_FFS = ["receptor_OpenMM_Gas", "receptor_OpenMM_OBC2", "Theta_1OpenMM_Gas", "Theta_1OpenMM_OBC2",
               "Theta_RLOpenMM_Gas", "Theta_RLOpenMM_OBC2", "Theta_1sander_Gas", "Theta_1sander_PBSA",
               "Theta_RLsander_Gas", "Theta_RLsander_PBSA", "receptor_sander_Gas", "receptor_sander_PBSA",
               "grid_MBAR", "MBAR", "OpenMM_Gas_MBAR_c2", "OpenMM_Gas_inverse_FEP",
               "OpenMM_OBC2_inverse_FEP", "OpenMM_OBC2_MBAR_c2"]


class MultiStruScores:
    """
    load and averaging scores for a ligand over multiple receptor's snapshots
    """
    def __init__( self, score_dir, ligand_group, ligand_3l_code, weights, yank_systems,
                yank_interaction_energies,
                exclude_ffs=EXCLUDE_FFS, 
                repeats=100):
        """
        :param score_dir: str
        :param ligand_group: str
        :param ligand_3l_code: str
        :param weights: dict, weights[system][snapshot] -> float
        :param yank_systems: list of str
        :param yank_interaction_energies:
        :param exclude_ffs: list of str
        :param repeats: int, number of bootstrap repeats
        """
        self._identification = ligand_group + ligand_3l_code
        self._dock_dir = os.path.join(score_dir, DOCK6_SUB_DIR, ligand_group, ligand_3l_code)
        self._algdock_dir = os.path.join(score_dir, ALGDOCK_SUB_DIR, ligand_group, ligand_3l_code)

        self._yank_interaction_energies = yank_interaction_energies
        
        self._FFs = ["dock6"]
        FFs = glob.glob(os.path.join(self._algdock_dir, "*.score"))
        FFs = [os.path.basename(FF)[:-6] for FF in FFs]
        FFs = [ff for ff in FFs if ff not in exclude_ffs]
        self._FFs.extend(FFs)
        
        self._scores = {}
        self._load_dock6()
        self._load_algdock()
         
        self._weights = weights
        self._yank_systems = yank_systems
        self._repeats = repeats
        self._considered_snapshots = self._get_snapshots()
        self._allowed_snapshots = self._get_snapshots_in_scores_and_systems()

    def _load_dock6(self):
        """
        store dock6 scores to self._scores["dock6"]
        self._scores["dock6"][snapshot_id] -> float
        :return: None
        """
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
        """
        store algdock scores into self._scores[FF], where FF != "dock6"
        self._scores[FF][snapshot_id] -> float
        :return: None
        """
        algdock_FFs = [FF for FF in self._FFs if FF != "dock6"]

        for FF in algdock_FFs:
            in_file = open(os.path.join(self._algdock_dir, FF + ".score"), "r")
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
        """
        :return: list of str
                list of all snapshots in  self._weights
        """
        snapshots = []
        for system in self._yank_systems:
            for snapshot in self._weights[system].keys():
                snapshots.append(snapshot)
        return snapshots

    def _get_snapshots_in_scores_and_systems(self):
        """
        :return:  snapshots, dict, snapshots[FF][system] -> list
                   snapshots in yank_systems and have scores
        """
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
        """
        If pass yank_systems with a single system, then we get relative binding free energy with to that system
        :param snapshots: lis of str
        :return: averages, dict, averages[FF] -> float
                -1/\beta *\ln(<e^{-\beta * BPMF}>),
                where <...> is average over snapshots
        """
        averages = {}
        for FF in self._FFs:
            sys_mean = 0.
            for system in self._yank_systems:
                a = 0.
                w = 0.
                for snapshot in snapshots:
                    if snapshot in self._allowed_snapshots[FF][system]:
                        try:
                            a += np.exp(-1.0 * (self._scores[FF][snapshot] -
                                                self._yank_interaction_energies[system][snapshot])) * self._weights[system][snapshot]
                        except FloatingPointError:
                            print("overflow for", self._identification, snapshot, FF)
                        else:
                            w += self._weights[system][snapshot]
                if w != 0:
                    a = a / w
                sys_mean += a * self._weights["systems"][system]

            sys_mean = sys_mean / sum([self._weights["systems"][system] for system in self._yank_systems])
            averages[FF] = (-1./BETA) * np.log(sys_mean)
        return averages

    def get_exp_mean(self):
        return self._cal_exp_mean(self._considered_snapshots)

    def cal_exp_mean_separate_for_each_system(self, FF, snapshots):
        """
        exponential mean for each FF, with respect to each system in self._yank_systems
        :param FF: str
        :param snapshots: list of str
        :return: sys_means, dict, sys_means[system] -> float
        """
        sys_means = {}
        for system in self._yank_systems:
            a = 0.
            w = 0.
            for snapshot in snapshots:
                if snapshot in self._allowed_snapshots[FF][system]:
                    try:
                        a += np.exp(-1.0 * (self._scores[FF][snapshot] -
                                            self._yank_interaction_energies[system][snapshot])) * self._weights[system][snapshot]
                    except FloatingPointError:
                        #print("overflow for", self._identification, snapshot, FF)
                        pass
                    else:
                        w += self._weights[system][snapshot]
            if w != 0:
                a = a / w
            a = (-1./BETA) * np.log(a)
            sys_means[system] = a

        return sys_means

    def get_exp_mean_std(self):
        """
        use bootstrap to estimate standard error of self._cal_exp_mean
        :return: std, dict
                 std[FF] -> float
        """
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
        """
        :param snapshots: list of str
        :return: averages, dict
                 averages[FF] -> float
        """
        averages = {}
        for FF in self._FFs:
            sys_mean = 0.
            for system in self._yank_systems:
                a = 0.
                w = 0.
                for snapshot in snapshots:
                    if snapshot in self._allowed_snapshots[FF][system]:
                        if self._scores[FF][snapshot] != np.inf:
                            a += (self._scores[FF][snapshot]
                                - self._yank_interaction_energies[system][snapshot])* self._weights[system][snapshot]
                            w += self._weights[system][snapshot]
                if w != 0:
                    a = a / w
                sys_mean += a * self._weights["systems"][system]

            sys_mean = sys_mean / np.sum([self._weights["systems"][system] for system in self._yank_systems])
            averages[FF] = sys_mean * TEMPERATURE * KB

        return averages

    def get_mean(self):
        return self._cal_mean(self._considered_snapshots)

    def get_mean_std(self):
        """
        use bootstrap to estimate standard error for self._cal_mean()
        :return: std, dict, std[FF] -> float
        """
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
        """
        :param snapshots: list of str
        :return: averages, dict
                 averages[FF] -> float
        """
        averages = {}
        for FF in self._FFs:
            a = []
            for system in self._yank_systems:
                a += [(self._scores[FF][snapshot] - self._yank_interaction_energies[system][snapshot])
                      for snapshot in snapshots if snapshot in self._allowed_snapshots[FF][system]]

            if len(a) > 0:
                averages[FF] = np.array(a).min() * TEMPERATURE * KB
            else:
                averages[FF] = np.inf

        return averages

    def get_min(self):
        return self._cal_min(self._considered_snapshots)

    def get_min_std(self):
        """
        use bootstrap for self._cal_min()
        :return: std, dict, std[FF] -> float
        """
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

    def check_extreme_low(self, cutoff=-100.):
        for FF in self._FFs:
            for system in self._yank_systems:
                for snapshot in self._weights[system].keys():
                    if snapshot in self._allowed_snapshots[FF][system]:
                        if self._scores[FF][snapshot] != np.inf:
                            if self._scores[FF][snapshot] < cutoff:
                                print("Extreme low: " + self._identification,
                                      snapshot, FF, " %20.10f"%self._scores[FF][snapshot])
        return None


