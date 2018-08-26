
import os
import argparse

from _algdock import SIX_YANK_SYSTEMS_NAMES

import sys
sys.path.append("/home/tnguye46/FFT_T4/scripts")
from _plots import scatter_plot, scatter_plot_info
from _yank_algdock_fft_scores import load_scores, matching_scores, write_pairs

parser = argparse.ArgumentParser()

parser.add_argument( "--score_dir",             type=str, default="/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/Relative_FE_vs_Yank_Est_1/subtract_self_fe/exclude_all_ref_ligands")
parser.add_argument( "--weight_scheme",         type=str, default="__equal_sys__single_weight")
parser.add_argument( "--score_file",            type=str, default="raw_data.dat")

parser.add_argument( "--out_figure",         type=str, default = "scatter_plot.pdf")
parser.add_argument( "--out_text",           type=str, default = "scatter_plot.dat")

args = parser.parse_args()

print "score_dir", args.score_dir
print "weight_scheme", args.weight_scheme

ref_ligands_code = SIX_YANK_SYSTEMS_NAMES.keys()
n_ref_ligands = len(ref_ligands_code)

for row in range(n_ref_ligands - 1):
    row_ligand_code = ref_ligands_code[row]
    row_ligand_name = SIX_YANK_SYSTEMS_NAMES[row_ligand_code]

    row_score_file = os.path.join(args.score_dir, row_ligand_code + args.weight_scheme, args.score_file)
    row_scores, row_stds = load_scores(row_score_file, 0, 2, 4, [])

    for col in range(row+1, n_ref_ligands):
        col_ligand_code = ref_ligands_code[col]
        col_ligand_name = SIX_YANK_SYSTEMS_NAMES[col_ligand_code]

        col_score_file = os.path.join(args.score_dir, col_ligand_code + args.weight_scheme, args.score_file)
        print "between " + row_ligand_name + " and " + col_ligand_name

        col_scores, col_stds = load_scores(col_score_file, 0, 2, 4, [])

        out_dir = row_ligand_name + "_vs_" + col_ligand_name
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        out_figure_file    = os.path.join(out_dir, args.out_figure)
        out_text_file      = os.path.join(out_dir, args.out_text)
        out_raw_data_file  = os.path.join(out_dir, "raw_data.dat")

        
        x, y, xerr, yerr, ligands = matching_scores(col_scores, row_scores, col_stds, row_stds)
        write_pairs(col_scores, row_scores, col_stds, row_stds, out_raw_data_file, [])

        # use one std for errorbar
        xerr /= 2.
        yerr /= 2.

        #markercolors = ["b" if ".inactive." in ligand or ligand == "phenol.A__AAA" else "r" for ligand in ligands]
        markers = ["D" if ".inactive." in ligand or ligand == "phenol.A__AAA" else "." for ligand in ligands]
        markercolors = ["k" for ligand in ligands]

        xlabel = col_ligand_name
        ylabel = row_ligand_name

        if len(x) > 1:
            scatter_plot_info(x, y, ligands, out_text_file)

            scatter_plot(x, y, xlabel, ylabel, out_figure_file, 
                        show_xy_axes=True,
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
                        text_pos=[0.1, 0.7] 
                        )

print "DONE"
