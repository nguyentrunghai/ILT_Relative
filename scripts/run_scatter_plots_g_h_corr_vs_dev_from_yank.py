"""
"""

import argparse

from _yank import load_scores

parser = argparse.ArgumentParser()

parser.add_argument("--yank_results", type=str, default="/home/tnguye46/T4_Lysozyme/Yank/yank_results.dat")

args = parser.parse_args()

yank_scores, yank_stds = load_scores(args.yank_results, 0, 1, 2, [])