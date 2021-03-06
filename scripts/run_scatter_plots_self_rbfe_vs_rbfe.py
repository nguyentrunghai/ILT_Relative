"""
run scatter plots between self relative binding free energies and relative binding free energies
TODO set figure titles
"""
from __future__ import print_function

import os
import glob
import argparse

import numpy as np

from _plots import scatter_plot

parser = argparse.ArgumentParser()

parser.add_argument("--data_dir", type=str,
default="Relative_FE_Est_with_CV_method_3a/flip_sign_c__not_subtract_self/benzene.A__AAA__equal_sys__single_weight/ExpMean")

parser.add_argument("--glob_matching", type=str, default="*_VERSUS_*")

parser.add_argument("--xlabel", type=str, default="Self relative binding free energy (kcal/mol)")
parser.add_argument("--ylabel", type=str, default="Relative binding free energy (kcal/mol)")

parser.add_argument("--remove_outliers", action="store_true", default=False)
parser.add_argument("--how_many_std", type=float, default=3.)

parser.add_argument("--title", action="store_true", default=False)

args = parser.parse_args()


def _iqr(x):
    q75, q25 = np.percentile(x, [75, 25])
    return q75 - q25


def _std_from_iqr(x):
    return _iqr(x) / 1.35


def _outliers(x, how_many_std=3):
    """
    :param x: 1d array
    :param how_many_std: float
    :return outliers: 1d bool array
    """
    x = np.asarray(x)
    assert x.ndim == 1, "x must be 1d"

    std = _std_from_iqr(x)
    median = np.median(x)

    lower = median - how_many_std * std
    upper = median + how_many_std * std

    outliers = (x < lower) | (x > upper)
    return outliers


def _remove_outliers(x, y, how_many_std=3):
    """
    :param x: 1d array
    :param y: 1d array
    :param weights: 1d array
    :return (new_x, new_y): 1d arrays, x, y after remove outliers in both
    """
    x = np.asarray(x)
    y = np.asarray(y)

    assert x.shape == y.shape, "x, y must have the same shape"

    outliers_x = _outliers(x, how_many_std=how_many_std)
    outliers_y = _outliers(y, how_many_std=how_many_std)
    all_outliers = outliers_x | outliers_y
    not_outliers = ~all_outliers
    return x[not_outliers], y[not_outliers]


data_files = glob.glob(os.path.join(args.data_dir, args.glob_matching))

for data_file in data_files:
    print("Processing file", data_file)
    data = np.loadtxt(data_file)
    x = data[:, 0]
    y = data[:, 1]

    if args.remove_outliers:
        print("Remove outliers beyond %0.2f std" % args.how_many_std)
        x, y = _remove_outliers(x, y, how_many_std=args.how_many_std)

    out_file = os.path.basename(data_file) + ".pdf"

    title = os.path.basename(data_file)
    if args.title:
        title = " vs ".join(title.split("_vs_"))
    else:
        title = None
    scatter_plot(x, y, args.xlabel, args.ylabel, out_file,
                 show_xy_axes=True,
                 xerr=None, yerr=None,
                 x_logscale=False,
                 y_logscale=False,
                 show_regression_line=False,
                 show_diagonal_line=False,
                 show_rmse=True,
                 show_R=True,
                 show_regression_line_eq=True,
                 markers=None,
                 markersize=4,
                 markercolors=None,
                 same_xy_scale=False,
                 integer_limits=False,
                 text_pos=[0.1, 0.7],
                 title=title)

