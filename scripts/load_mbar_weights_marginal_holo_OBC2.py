
import os
import numpy as np

from _algdock import SIX_YANK_SYSTEMS

def load_mbar_weights(path="/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/Mbar_Weights_marginal_holo_OBC2"):
    """
    :param path: str
    :return: weights, dict, weights[ligand][snapshot] -> float
    """
    weights = dict()
    weights["systems"] = {ligand:1 for ligand in SIX_YANK_SYSTEMS}

    for ligand in SIX_YANK_SYSTEMS:
        in_file = os.path.join(path, ligand + ".weig")
        weights[ligand] = {}
        with open(in_file, "r") as handle:
            for line in handle:
                snapshot, weight = line.split()
                weights[ligand][snapshot] = np.float(weight)
    return weights

