
from __future__ import print_function

import argparse
import glob
import os
import numpy as np

from _algdock import SIX_YANK_SYSTEMS

import sys
sys.path.append("/home/tnguye46/FFT_T4/scripts")
from _plots import scatter_plot, scatter_plot_info
from _yank_algdock_fft_scores import load_scores, matching_scores, write_pairs


parser = argparse.ArgumentParser()
parser.add_argument("--yank_results", type=str, default="/home/tnguye46/T4_Lysozyme/Yank/yank_results.dat")
parser.add_argument("--algdock_results_dir", type=str,
                    default="/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/Relative_FE_Est_1/OBC2")
parser.add_argument("--averaging_rule", type=str, default="ExpMean")
parser.add_argument("--algdock_score_file", type=str, default="OpenMM_OBC2_MBAR.score")

parser.add_argument("--exclude_ligands_from_scatter_plots", type=str, default="")

parser.add_argument("--subtract_self_rbfe", action="store_true", default=False)
parser.add_argument("--compare_absolute",  action="store_true", default=False)

parser.add_argument("--out_figure", type=str, default="algdock_vs_yank.pdf")
parser.add_argument("--out_text", type=str, default="algdock_vs_yank.dat")

parser.add_argument("--xlabel", type=str, default="YANK free energy (kcal/mol)")
parser.add_argument("--ylabel", type=str, default="AlGDock free energy (kcal/mol)")
parser.add_argument("--show_xy_axes", type=bool, default=True)

args = parser.parse_args()

exclude_ligands_from_scatter_plots = args.exclude_ligands_from_scatter_plots.split()
print("exclude_ligands_from_scatter_plots ", exclude_ligands_from_scatter_plots)
print("subtract_self_rbfe ", args.subtract_self_rbfe)
print("compare_absolute ", args.compare_absolute)

yank_scores, yank_stds = load_scores(args.yank_results, 0, 1, 2, [])

weighting_schemes_ref_ligands = glob.glob(os.path.join(args.algdock_results_dir, "*",
                                                       args.averaging_rule, args.algdock_score_file))
weighting_schemes_ref_ligands = [file_name.split("/")[-3] for file_name in weighting_schemes_ref_ligands]

weighting_schemes = []
for scheme_ligand in weighting_schemes_ref_ligands:
    for ligand in SIX_YANK_SYSTEMS:
        if ligand in scheme_ligand:
            scheme = scheme_ligand.split(ligand)[-1]
            weighting_schemes.append(scheme)

weighting_schemes = list(set(weighting_schemes))
print(weighting_schemes)

all_algdock_scores = {}
all_algdock_stds = {}
for scheme in weighting_schemes:
    all_algdock_scores[scheme] = {}
    all_algdock_stds[scheme] = {}

    for ref_ligand in SIX_YANK_SYSTEMS:

        out_dir = ref_ligand + scheme
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        out_figure_file = os.path.join(out_dir, args.out_figure)
        out_text_file = os.path.join(out_dir, args.out_text)
        out_raw_data_file = os.path.join(out_dir, "raw_data.dat")

        algdock_score_file = os.path.join(args.algdock_results_dir, out_dir,
                                          args.averaging_rule, args.algdock_score_file)
        print(algdock_score_file)
        algdock_scores, algdock_stds = load_scores(algdock_score_file, 0, 1, 2, [])

        if args.subtract_self_rbfe:
            algdock_scores = {ligand: (algdock_scores[ligand] - algdock_scores[ref_ligand])
                              for ligand in algdock_scores}

        if args.compare_absolute:
            algdock_scores = {ligand: (algdock_scores[ligand] + yank_scores[ref_ligand]) for ligand in algdock_scores}

        
        all_algdock_scores[scheme][ref_ligand] = algdock_scores
        all_algdock_stds[scheme][ref_ligand] = algdock_stds


        considered_algdock_scores = {ligand: algdock_scores[ligand] for ligand in algdock_scores
                                     if ligand not in [ref_ligand] + exclude_ligands_from_scatter_plots}

        if args.compare_absolute:
            considered_yank_scores = yank_scores
        else:
            considered_yank_scores = {ligand: (yank_scores[ligand] - yank_scores[ref_ligand]) for ligand in yank_scores}

        x, y, xerr, yerr, ligands = matching_scores(considered_yank_scores, considered_algdock_scores,
                                                    yank_stds, algdock_stds)
        write_pairs(considered_yank_scores, considered_algdock_scores, yank_stds, algdock_stds, out_raw_data_file, [])
        # use one std for errorbar
        xerr /= 2.
        yerr /= 2.

        #markercolors = ["b" if ".inactive." in ligand or ligand == "phenol.A__AAA" else "r" for ligand in ligands]
        markers = ["D" if ".inactive." in ligand or ligand == "phenol.A__AAA" else "." for ligand in ligands]
        markercolors = ["k" for ligand in ligands]

        if len(x) > 1:
            scatter_plot_info(x, y, ligands, out_text_file)

            scatter_plot(x, y, args.xlabel, args.ylabel, out_figure_file, 
                        show_xy_axes=args.show_xy_axes,
                        xerr=xerr, yerr=yerr,
                        show_regression_line=True,
                        show_diagonal_line=False,
                        show_rmse=True,
                        show_R=True,
                        show_regression_line_eq=True,
                        markers=markers,
                        markersize=5,
                        markercolors=markercolors,
                        same_xy_scale=False,
                        text_pos=[0.1, 0.7])


all_ligands = [all_algdock_scores[scheme][ref_ligand].keys()
               for scheme in weighting_schemes for ref_ligand in SIX_YANK_SYSTEMS]
common_ligands = set(all_ligands[0])
unioned_ligands = set(all_ligands[0])

for ligand_set in all_ligands:
    common_ligands = common_ligands.intersection(ligand_set)
    unioned_ligands = unioned_ligands.union(ligand_set)

common_ligands = list(common_ligands)

min_algdock_scores = {}
min_algdock_stds = {}
for scheme in weighting_schemes:
    min_algdock_scores[scheme] = {}
    min_algdock_stds[scheme] = {}
    for ligand in common_ligands:
        min_algdock_scores[scheme][ligand] = np.min([all_algdock_scores[scheme][ref_ligand][ligand]
                                                     for ref_ligand in SIX_YANK_SYSTEMS])
        min_algdock_stds[scheme][ligand] = np.mean([all_algdock_stds[scheme][ref_ligand][ligand]
                                                    for ref_ligand in SIX_YANK_SYSTEMS])

# take all 24 ligands

for scheme in weighting_schemes:
    min_algdock_scores[scheme] = {}
    min_algdock_stds[scheme] = {}

    for ligand in unioned_ligands:
        m_algdock_score = []
        m_algdock_std = []
        for ref_ligand in SIX_YANK_SYSTEMS:
            if ref_ligand in all_algdock_scores[scheme]:
                if ligand in all_algdock_scores[scheme][ref_ligand]:

                    m_algdock_score.append( all_algdock_scores[scheme][ref_ligand][ligand] )
                    m_algdock_std.append( all_algdock_stds[scheme][ref_ligand][ligand] )
        #
        min_algdock_scores[scheme][ligand] = np.min(m_algdock_score)
        min_algdock_stds[scheme][ligand] = np.mean(m_algdock_std)


for scheme in weighting_schemes:
    out_dir = "min" + scheme
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    out_figure_file = os.path.join(out_dir, args.out_figure)
    out_text_file = os.path.join(out_dir, args.out_text)
    out_raw_data_file = os.path.join(out_dir, "raw_data.dat")

    if args.compare_absolute:
        considered_yank_scores = yank_scores
    else:
        considered_yank_scores = {ligand: (yank_scores[ligand] - yank_scores[ref_ligand]) for ligand in yank_scores}

    considered_yank_scores = {ligand: considered_yank_scores[ligand] for ligand in considered_yank_scores
                              if ligand not in exclude_ligands_from_scatter_plots}

    x, y, xerr, yerr, ligands = matching_scores(considered_yank_scores, min_algdock_scores[scheme],
                                                yank_stds, min_algdock_stds[scheme])
    write_pairs(considered_yank_scores, min_algdock_scores[scheme], yank_stds,
                min_algdock_stds[scheme], out_raw_data_file, [])
    # use one std for errorbar
    xerr /= 2.
    yerr /= 2.

    #markercolors = ["b" if ".inactive." in ligand or ligand == "phenol.A__AAA" else "r" for ligand in ligands]
    markers = ["D" if ".inactive." in ligand or ligand == "phenol.A__AAA" else "." for ligand in ligands]
    markercolors = ["k" for ligand in ligands]

    if len(x) > 1:
        scatter_plot_info(x, y, ligands, out_text_file)

        scatter_plot(x, y, args.xlabel, args.ylabel, out_figure_file, 
                    show_xy_axes=args.show_xy_axes,
                    xerr=xerr, yerr=yerr,
                    show_regression_line=True,
                    show_diagonal_line=False,
                    show_rmse=True,
                    show_R=True,
                    show_regression_line_eq=True,
                    markers=markers,
                    markersize=5,
                    markercolors=markercolors,
                    same_xy_scale=False,
                    text_pos=[0.1, 0.7])


print("Done")

