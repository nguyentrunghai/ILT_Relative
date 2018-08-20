

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import seaborn
seaborn.set_style("white")

from _algdock import SIX_YANK_SYSTEMS


def _transform_weight(weights, max_to=100):
    w_max = weights.max()
    d = max_to / w_max
    return weights * d

def _plot(hx, hy, x, y, weights,
        xlabel, ylabel, out,
        xlimits=None,
        ylimits=None,
        nticks=8,
        figure_size=(3.2, 3.2*6/8),
        dpi=300,
        fontsize=8,
        font = {"fontname": "Arial"}):
    """
    """
    my_cmap = plt.get_cmap('Reds')
    plt.figure(figsize=figure_size)
    plt.hist2d(hx, hy, bins=20, cmap=my_cmap, normed=True)
    cbar = plt.colorbar(fraction=0.046, pad=0.04)

    plt.scatter(x, y, s=weights, alpha=0.6, facecolors='none', edgecolors='k')

    ax = plt.axes()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)

    plt.xlabel(xlabel, fontsize=fontsize, **font)
    plt.ylabel(ylabel, fontsize=fontsize, **font)

    cbar.locator = ticker.MaxNLocator(nbins=nticks)
    cbar.update_ticks()
    cbar.ax.tick_params(labelsize=fontsize)

    ax.locator_params(axis='x', nbins=nticks)
    ax.locator_params(axis='y', nbins=nticks)

    if xlimits is not None:
        ax.set_xlim(xlimits) 

    if ylimits is not None:
        ax.set_ylim(ylimits) 

    plt.tight_layout()
    plt.savefig(out, dpi=dpi)
    return None


pc_files = {yank_system : yank_system + "_pc.dat" for yank_system in SIX_YANK_SYSTEMS}
pc_files["all_24_holo"] = "all_24_holo.dat"

pc_data = { key : np.loadtxt( pc_files[key] ) for key in pc_files }

x_min = np.min( [ pc_data[key][:,0].min() for key in pc_data] )
x_max = np.max( [ pc_data[key][:,0].max() for key in pc_data] )

y_min = np.min( [ pc_data[key][:,1].min() for key in pc_data] )
y_max = np.max( [ pc_data[key][:,1].max() for key in pc_data] )

xlabel = "PC1"                                                                                                         
ylabel = "PC2"

hx = pc_data["all_24_holo"][:,0]
hy = pc_data["all_24_holo"][:,1]

for ligand in SIX_YANK_SYSTEMS:
    x = pc_data[ligand][:, 0]
    y = pc_data[ligand][:, 1]
    weights = pc_data[ligand][:, 2]
    weights = _transform_weight(weights)
    out = ligand + ".pdf"
    _plot(hx, hy, x, y, weights, xlabel, ylabel, out,
            xlimits=[x_min, x_max],
            ylimits=[y_min, y_max])

