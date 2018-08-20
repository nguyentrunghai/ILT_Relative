
import sys
import pickle
import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import seaborn as sns
sns.set_style("whitegrid")

from _yank import YANK_LIGANDS

sys.path.append("/home/tnguye46/FFT_T4/scripts")
from _plots import scatter_plot, scatter_plot_info
from _yank_algdock_fft_scores import load_scores

REF_HOLO_ENSEMBLES = ["1-methylpyrrole.A__AAA", "benzene.A__AAA", "p-xylene.A__AAA",
                            "phenol.A__AAA", "lysozyme.inactive.A__AAS", "lysozyme.active.A__ABJ"]
ALL_YANK_LIGANDS = YANK_LIGANDS.keys()
TARGET_LIGANDS_18 = [ligand for ligand in ALL_YANK_LIGANDS if ligand not in REF_HOLO_ENSEMBLES]

TEMPERATURE = 300.0
KB = 8.3144621E-3/4.184  # kcal/mol/K
BETA = 1. / TEMPERATURE / KB
V_0 = 1661.

def algdock_site_correction(r0):
    V_binding = 4. / 3. * np.pi * (r0 ** 3)
    return -TEMPERATURE * KB * np.log(V_binding / V_0 / 8 / np.pi**2)

SITE_CORRECTION = algdock_site_correction( 5.0 )


def _load_yank_abs_fe(yank_result_file):
    yank_abs_fe, _ = load_scores(yank_result_file, 0, 1, 2, [])
    return yank_abs_fe


def _load_algdock(algdock_score_file):
    algdock_score, _ = load_scores(algdock_score_file, 0, 1, 2, [])
    return algdock_score


def _yank_relative_fe_matrix(yank_abs_fe):
    """
    return a matrix, where col are 6 ref holo ensembles, rows are 18 target holo ensembles
    """
    yank_relat_fe_matrix = {}

    for ref_ligand in REF_HOLO_ENSEMBLES:
        yank_relat_fe_matrix[ref_ligand] = {}

        for target_ligand in TARGET_LIGANDS_18:
            yank_relat_fe_matrix[ref_ligand][target_ligand] = yank_abs_fe[target_ligand] - yank_abs_fe[ref_ligand]
    return yank_relat_fe_matrix


def _load_algdock_relat_fe_matrix(algdock_relat_fe_files):
    """
    algdock_relat_fe_files  : dict mapping reference ligand to algdock relative fe file
    return a matrix, where col are 6 ref holo ensembles, rows are 18 target holo ensembles
    """
    algdock_relat_fe_matrix = {}

    for ref_ligand in REF_HOLO_ENSEMBLES:
        algdock_relat_fe_matrix[ref_ligand] = {}
        relat_fes = _load_algdock( algdock_relat_fe_files[ref_ligand] )

        for target_ligand in TARGET_LIGANDS_18:
            # subtract self fe
            algdock_relat_fe_matrix[ref_ligand][target_ligand] = relat_fes[target_ligand] - relat_fes[ref_ligand]
    return algdock_relat_fe_matrix


def _load_algdock_abs_fes(algdock_abso_fe_file):
    """
    """
    algdock_abs_fes = _load_algdock(algdock_abso_fe_file)
    for ligand in algdock_abs_fes:
        algdock_abs_fes[ligand] += SITE_CORRECTION
    return algdock_abs_fes


def _load_jsd_matrix(jsd_matrix_file):
    return pickle.load(open(jsd_matrix_file, "r")) 



def _scatter_plot(xs, ys, out, 
                    xlimit=None,
                    ylimit=None,
                    xlabel=None, ylabel=None,
                    label_fontsize=8,
                    markers=None,
                    markersize=20,
                    facecolors=None, 
                    edgecolors='none',
                    alpha=1,
                    legends=None, 
                    legend_pos="best", 
                    legend_fontsize=6,
                    ntics=10,
                    tick_fontsize=8,

                    figure_size=(3.2, 2.4),
                    dpi=300,
                    font = {"fontname": "Arial"}
                    ):
    """
    """
    fig = plt.figure(figsize=figure_size)
    ax = plt.axes()

    if facecolors is None:
        facecolors = [ None for _ in range(len(xs)) ]

    if edgecolors is None:
        edgecolors = [ 'none' for _ in range(len(xs)) ]
    elif edgecolors == "face":
        edgecolors = [ 'face' for _ in range(len(xs)) ]
    else:
        raise Exception("unknown edgecolors")

    if markers is None:
        markers = [ "o" for _ in range(len(xs)) ]

    x_min = np.min([np.min(x) for x in xs])
    x_max = np.max([np.max(x) for x in xs])
    #ax.plot([x_min, x_max], [0, 0], color="k", linestyle='-', lw=1, alpha=alpha)

    scat_plots = []
    for i in range(len(xs)):
        pl = ax.scatter(xs[i], ys[i], s=markersize, marker=markers[i], alpha=alpha, facecolors=facecolors[i], edgecolors=edgecolors[i])
        scat_plots.append(pl)

    if xlimit is not None:
        ax.set_xlim(xlimit)
        
    if ylimit is not None:
        ax.set_ylim(ylimit)

    ax.locator_params(axis='x', nbins=ntics)
    ax.locator_params(axis='y', nbins=ntics)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(tick_fontsize)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(tick_fontsize)

    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=label_fontsize, **font)

    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=label_fontsize, **font)

    if legends is not None:
        ax.legend(scat_plots, legends , loc=legend_pos, fancybox=False, fontsize=legend_fontsize)
        #ax.legend(legends, loc='lower left', fancybox=True, ncol=3, fontsize=legend_fontsize)

    plt.tight_layout()
    plt.savefig(out, dpi=dpi)
    return None

parser = argparse.ArgumentParser()
parser.add_argument( "--yank_results",  type=str, default="/home/tnguye46/T4_Lysozyme/Yank/yank_results.dat")

parser.add_argument( "--algdock_relat_fe_dir",  type=str, default="/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/Relative_FE_Est_1/6holo")
parser.add_argument( "--algdock_abs_fe_dir",    type=str, default="/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/Absolute_FE/6apo/all_snap__equal_sys__state_weight/ExpMean")

parser.add_argument( "--algdock_score_file",         type=str, default = "OpenMM_OBC2_MBAR.score")

parser.add_argument( "--jsd_matrix_file",    type=str, default="/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/pca_apo_vs_holo_density/b0.1/jsd_matrix_ligand_code.pkl")

parser.add_argument( "--xlimit",    type=str, default=[0, 0.6])
parser.add_argument( "--ylimit",    type=str, default=[-6, 10])

parser.add_argument( "--out",    type=str, default="dev_from_yank_vs_jsd_6holo_6apo.pdf")
args = parser.parse_args()


colors = {"apo":"r"}
colors.update( { ref_ligand : c for ref_ligand, c in zip(REF_HOLO_ENSEMBLES, ("b", "g", "c", "m", "y", "k")) } )

markers = {"apo":"o"}
markers.update( { ref_ligand : c for ref_ligand, c in zip(REF_HOLO_ENSEMBLES, ("v", "^", "<", ">", "s", "D")) } )

legends = {"apo":"apo"}
legends.update( { ref_ligand:YANK_LIGANDS[ref_ligand] for ref_ligand in REF_HOLO_ENSEMBLES } )

yank_abs_fes = _load_yank_abs_fe(args.yank_results)
yank_relat_fe_matrix = _yank_relative_fe_matrix(yank_abs_fes)

algdock_relat_fe_files = { ref_ligand : os.path.join(args.algdock_relat_fe_dir, ref_ligand+"__equal_sys__single_weight", "ExpMean", args.algdock_score_file) 
                        for ref_ligand in REF_HOLO_ENSEMBLES}
algdock_relat_fe_matrix = _load_algdock_relat_fe_matrix(algdock_relat_fe_files)

algdock_abs_fes = _load_algdock_abs_fes( os.path.join( args.algdock_abs_fe_dir, args.algdock_score_file ) )

jsd_matrix = _load_jsd_matrix(args.jsd_matrix_file) 


xs = []
ys = []
cs = []
ms = []
ls = []

xs.append( [ jsd_matrix["apo"][target_ligand] for target_ligand in TARGET_LIGANDS_18 ] )
ys.append( [ algdock_abs_fes[target_ligand] - yank_abs_fes[target_ligand] for target_ligand in TARGET_LIGANDS_18 ] )
cs.append(colors["apo"])
ms.append(markers["apo"])
ls.append(legends["apo"])


for ref_ligand in REF_HOLO_ENSEMBLES:

    xs.append( [ jsd_matrix[ref_ligand][target_ligand] for target_ligand in TARGET_LIGANDS_18 ] )

    ys.append( [ algdock_relat_fe_matrix[ref_ligand][target_ligand] - yank_relat_fe_matrix[ref_ligand][target_ligand] for target_ligand in TARGET_LIGANDS_18 ] )
    
    cs.append(colors[ref_ligand])
    
    ms.append(markers[ref_ligand])

    ls.append(legends[ref_ligand])

_scatter_plot(xs, ys, args.out,
                xlimit=args.xlimit,
                ylimit=args.ylimit,
                xlabel="$D_{JS}$", 
                ylabel="Deviation from YANK (kcal/mol)",
                markers=ms,
                facecolors=cs, 
                edgecolors='face',
                alpha=0.5,
                legends=None
                ) 

print "DONE"
