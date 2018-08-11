
import argparse
import numpy as np

import sys
sys.path.append("/home/tnguye46/FFT_T4/scripts")
from _plots import improved_plot_lines

parser = argparse.ArgumentParser()
parser.add_argument( "--absolute_result",       type=str, default = "convergence_absolute_average_across_systems_R.dat")
#parser.add_argument( "--absolute_result",       type=str, default = "convergence_absolute_average_across_systems_RMSE.dat")

parser.add_argument( "--relative_result",       type=str, default = "convergence_relative_est_1_R.dat")
#parser.add_argument( "--relative_result",       type=str, default = "convergence_relative_est_1_RMSE.dat")

parser.add_argument( "--xlabel",       type=str, default = "# receptor snapshots")

#parser.add_argument( "--ylabel",       type=str, default = "Pearson's R w.r.t. YANK")
#parser.add_argument( "--ylabel",       type=str, default = "RMSE w.r.t. YANK (kcal/mol)")

parser.add_argument( "--ylabel",       type=str, default = "Pearson's R w.r.t. final results")
#parser.add_argument( "--ylabel",       type=str, default = "RMSE w.r.t. final results (kcal/mol)")

parser.add_argument( "--out",       type=str, default = "convergence_min_rel_average_abs_R.pdf")
#parser.add_argument( "--out",       type=str, default = "convergence_min_rel_average_abs_RMSE.pdf")

#parser.add_argument( "--out",       type=str, default = "convergence_min_rel_average_abs_not_convert_2_abs_R.pdf")
#parser.add_argument( "--out",       type=str, default = "convergence_min_rel_average_abs_not_convert_2_abs_RMSE.pdf")

#parser.add_argument( "--legend_pos",       type=str, default = "lower right")
parser.add_argument( "--legend_pos",       type=str, default = "upper right")

args = parser.parse_args()

absolute_data = np.loadtxt(args.absolute_result)
relative_data = np.loadtxt(args.relative_result)

xs    = [ absolute_data[:,0], relative_data[:,0] ]
ys    = [ absolute_data[:,1], relative_data[:,1] ]
yerrs = [absolute_data[:,2]/2., relative_data[:,2]/2. ]

legends=["Abs", "Rel"]

colors = ["b", "g"]
line_styles = ["-", "--"]
lw = 2.

#improved_plot_lines(xs, ys, yerrs=yerrs, xlabel=args.xlabel, ylabel=args.ylabel, out=args.out,
#                    legends=["Abs", "Rel"], legend_pos=args.legend_pos, 
#                    integer_xtics=True)

# no legend
improved_plot_lines(xs, ys, yerrs=yerrs, xlabel=args.xlabel, ylabel=args.ylabel, out=args.out,
                    integer_xtics=True,
                    colors=colors,
                    line_styles=line_styles,
                    lw=lw)

