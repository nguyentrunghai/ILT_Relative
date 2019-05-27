"""
implement estimator for relative binding free energies using holo estimator
with control variates.
The self relative free energy is used as a control variable.
"""
import numpy as np

from _relative_estimators_without_cv import MultiStruScores


def bootstrap_estimate_cov_var(score_dir, target_ligand, ref_ligand, all_ref_ligands,
            weights, yank_interaction_energies,
            FF, sample_size, repeats):
    """
    calculate Cov(relative_bfe, self_relative_bfe) and Var(self_relative_bfe)
    using bootstrap

    :param score_dir: str
    :param target_ligand: str
    :param ref_ligand: str
    :param all_ref_ligands: list of str, all possible ref ligand
    :param weights: dict,
                    weights[ref_ligand_name][snapshot] -> float
                    weights["systems"][ref_ligand_name] -> float
    :param yank_interaction_energies: dict, yank_interaction_energies[system][snapshot] -> float
    :param FF: str, phase
    :param sample_size: int
    :param repeats: int
    :return: (cov, var)
            cov, float, covariance
            var, float, variance
    """
    assert ref_ligand in all_ref_ligands, "Unknown ref ligand: "+ ref_ligand

    snapshots_for_ref_ligand = weights[ref_ligand].keys()
    assert sample_size <= len(snapshots_for_ref_ligand), "sample_size is larger than number of snapshots"

    target_ligand_group = target_ligand[:-3]
    target_ligand_3l_code = target_ligand[-3:]

    ref_ligand_group = ref_ligand[: -3]
    ref_ligand_3l_code = ref_ligand[-3:]

    rel_fes = []
    self_rel_fes = []
    for _ in range(repeats):
        random_snapshots = np.random.choice(snapshots_for_ref_ligand, size=sample_size, replace=True)

        rel_fe = MultiStruScores(score_dir, target_ligand_group, target_ligand_3l_code,
                                 weights, [ref_ligand], yank_interaction_energies
                                 ).cal_exp_mean_separate_for_each_system(FF, random_snapshots)

        self_rel_fe = MultiStruScores(score_dir, ref_ligand_group, ref_ligand_3l_code,
                                      weights, [ref_ligand], yank_interaction_energies
                                      ).cal_exp_mean_separate_for_each_system(FF, random_snapshots)

        if (rel_fe not in [np.nan, np.inf, -np.inf]) and (self_rel_fe not in [np.nan, np.inf, -np.inf]):
            rel_fes.append(rel_fe)
            self_rel_fes.append(self_rel_fe)

    return np.cov(rel_fes, self_rel_fes)[0, -1], np.var(self_rel_fes)

