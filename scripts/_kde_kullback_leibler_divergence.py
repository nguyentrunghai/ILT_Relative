
import numpy as np
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
from matplotlib import ticker
import seaborn as sns


def _mins(pcs_file):
    """
    :param pcs_file: list of str
    :return: mins, numpy array
    """
    return np.loadtxt(pcs_file).min(axis=0)


def _maxs(pcs_file):
    """
    :param pcs_file: list of str
    :return: maxs, numpy array
    """
    return np.loadtxt(pcs_file).max(axis=0)


def overal_min_max(pcs_files):
    """
    :param pcs_files: list str
    :return: min_max, list of tuples
             eg., [(pc1_min, pc1_max), (pc2_min, pc2_max), ...]
    """
    assert isinstance(pcs_files, list), "pcs_files must be a list of str"
    mins = [_mins(pcs_file) for pcs_file in pcs_files]
    maxs = [_maxs(pcs_file) for pcs_file in pcs_files]

    mins = np.array(mins).min(axis=0)
    maxs = np.array(maxs).max(axis=0)

    # add extra buffer
    factor = 0.00001
    min_max = []
    for mi, ma in zip(mins, maxs):
        d = ma - mi
        min_max.append( (mi - factor*d, ma + factor*d) )

    return min_max


def gaussian_kde_2d(x, y, xmin, xmax, ymin, ymax, bandwidth, xbins=20, ybins=20):
    """ Build 2D kernel density estimate (KDE)."""

    # xy_train with shape (nsamples, nfeatures)
    xy_train = np.hstack( [ x[:, np.newaxis], y[:, np.newaxis] ] )

    kde_skl = KernelDensity(bandwidth=bandwidth, kernel='gaussian')
    kde_skl.fit(xy_train)

    # x moves along axis 1
    # y moves along axis 0
    x_grid, y_grid = np.meshgrid( np.linspace(xmin, xmax, xbins), np.linspace(ymin, ymax, ybins) )
    xy_sample = np.hstack( [ x_grid.ravel()[:, np.newaxis], y_grid.ravel()[:, np.newaxis] ] )

    density_grid = np.exp( kde_skl.score_samples(xy_sample) )
    density_grid = density_grid.reshape(x_grid.shape)

    return x_grid, y_grid, density_grid


def bin_area(x_grid, y_grid):
    """
    :param x_grid: numpy array
    :param y_grid: numpy array
    :return: float
    """
    x_bin_width = x_grid[0, 1] - x_grid[0, 0]
    y_bin_width = y_grid[1, 0] - y_grid[0, 0]
    return x_bin_width * y_bin_width 


def plot_density(density_grid, left, right, bottom, top, 
                xlabel, ylabel, out,
                nticks=7):
    """
    :param density_grid:
    :param left:
    :param right:
    :param bottom:
    :param top:
    :param xlabel:
    :param ylabel:
    :param out:
    :param nticks:
    :return:
    """
    figure_size = (3.2, 3.2*6/8)
    dpi = 300
    fontsize = 8
    markersize = 8
    font = { "fontname": "Arial"}

    plt.figure(figsize=figure_size)

    # for more colormaps see http://matplotlib.org/examples/color/colormaps_reference.html
    my_cmap = plt.get_cmap('hot')

    plt.imshow(density_grid, cmap=my_cmap, interpolation='nearest', origin='lower', 
                extent=(left, right, bottom, top))

    ax = plt.axes()
    cbar = plt.colorbar(fraction=0.046, pad=0.04)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)

    ax.locator_params(axis='x', nbins=nticks)
    ax.locator_params(axis='y', nbins=nticks)

    cbar.locator = ticker.MaxNLocator(nbins=nticks)
    cbar.update_ticks()
    cbar.ax.tick_params(labelsize=fontsize)

    plt.xlabel(xlabel, fontsize=fontsize, **font)
    plt.ylabel(ylabel, fontsize=fontsize, **font)

    plt.tight_layout()
    plt.savefig(out, dpi=dpi)
    return None


def kullback_leibler_divergence(p, q, bin_area):
    """
    :param p:
    :param q:
    :param bin_area:
    :return:
    """
    assert p.shape == q.shape, "p and q must have the same shape"
    assert p.ndim == 2, "p must be 2d array"

    nx, ny = p.shape
    total = 0
    
    for i in range(nx):
        for j in range(ny):
            pc = p[i, j]
            qc = q[i, j]

            if(pc > 0 and qc > 0): 
                total += pc * (np.log(pc) - np.log(qc))

            elif(pc !=0 and qc == 0):
                total += pc * np.log(pc)

    return total * bin_area


def jensen_shannon_divergence(p, q, bin_area):
    m = 0.5 * (p + q)
    jsd = 0.5 * kullback_leibler_divergence(p, m, bin_area) + 0.5 * kullback_leibler_divergence(q, m, bin_area)
    return jsd


def annotated_heat_map(dataframe, out):
    """
    :param dataframe:
    :param out:
    :return:
    """
    figure_size=(6.4, 6.4)
    dpi=300
    sns.set(font_scale=1)
    
    plt.figure(figsize=figure_size)
    my_cmap = plt.get_cmap('Reds')
    
    ax = plt.axes()
    sns.heatmap(dataframe, annot=True, fmt=".2f", linewidths=.5,  ax=ax, cmap=my_cmap)
    
    row_names = list(dataframe.index)
    col_names = list(dataframe.columns)

    ax.set_yticklabels(row_names, rotation=0)
    ax.set_xticklabels(col_names, rotation=45)

    plt.tight_layout()
    plt.savefig(out, dpi=dpi)
    return

