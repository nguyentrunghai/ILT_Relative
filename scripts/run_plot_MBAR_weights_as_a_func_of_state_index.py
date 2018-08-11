
import sys
import numpy as np
from pandas import Series

from _algdock import SIX_YANK_SYSTEMS, load_algdock_snapshots_for_each_of_six_yank_systems
from load_mbar_weights_holo_OBC2 import load_mbar_weights

sys.path.append("/home/tnguye46/FFT_T4/scripts")
from _plots import improved_plot_lines

def _group_mean(weights, ordered_snapshots):
    w = Series(weights)
    w = w.reindex(ordered_snapshots)

    nstates = 16
    snapshots_per_state = len(w)/nstates

    state_indices = []
    for i in range(nstates):
        state_indices += [i] * snapshots_per_state
    assert len(state_indices) == len(w)

    g_mean = w.groupby(state_indices).mean()
    states = []
    mean_w = []
    for i, v in g_mean.iteritems():
        states.append(i)
        mean_w.append(v)

    return np.array(states, dtype=int), np.array(mean_w, dtype=float)


# take single_snap_weights
block_weights, state_weights, single_snap_weights, stru_group_weights_equal_sys, stru_group_weights_ub_weighted = load_mbar_weights()

algdock_snapshots_for_each_of_six_yank_systems = load_algdock_snapshots_for_each_of_six_yank_systems()

ligand_names = {}
ligand_names["lysozyme.active.A__ABJ"] = "n-hexylbenzene"
ligand_names["benzene.A__AAA"] = "benzene"
ligand_names["lysozyme.inactive.A__AAS"] = "($\pm$)-camphor"
ligand_names["1-methylpyrrole.A__AAA"] = "methylpyrrole"
ligand_names["phenol.A__AAA"] = "phenol"
ligand_names["p-xylene.A__AAA"] = "p-xylene"

legends = [ligand_names[ligand] for ligand in SIX_YANK_SYSTEMS]

xs = []
ys = []
for ligand in SIX_YANK_SYSTEMS:
    print ligand
    weights = single_snap_weights[ligand]
    ordered_snapshots = algdock_snapshots_for_each_of_six_yank_systems[ligand]

    states, mean_w = _group_mean(weights, ordered_snapshots)

    xs.append(states+1)
    ys.append(mean_w)

xlabel = "Alchemical state index"
ylabel = "Mean holo weight"

out = "holo_mbar_weights_vs_state_index.pdf"
improved_plot_lines(xs, ys, xlabel=xlabel, ylabel=ylabel, out=out,
        legends=legends,
        y_logscale=False)

