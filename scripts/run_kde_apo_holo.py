
import argparse
import glob
import os
import pickle

import numpy as np

from _kde_kullback_leibler_divergence import overal_min_max, gaussian_kde_2d, plot_density

parser = argparse.ArgumentParser()

parser.add_argument( "--pc_file_glob_matching",  type=str, default="/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/pca_apo_vs_holo/*.dat")

parser.add_argument( "--bandwidth",  type=float, default=0.1)
parser.add_argument( "--xbins",  type=str, default=25)
parser.add_argument( "--ybins",  type=str, default=25)

parser.add_argument( "--out_prefix",  type=str, default="density")

args = parser.parse_args()

pc_files = glob.glob(args.pc_file_glob_matching)

(pc1_min, pc1_max), (pc2_min, pc2_max) = overal_min_max(pc_files)

for pc_file in pc_files:
    
    base_name = os.path.basename(pc_file)
    
    if "apo" in base_name:
        name = base_name
    elif "holo" in base_name:
        name = base_name.split("holo_")[-1]
    else:
        raise Exception("Unknown file: " + base_name)

    name = name.split(".dat")[0]
    print name 

    pc1_2 = np.loadtxt(pc_file)

    pc1_grid, pc2_grid, density_grid = gaussian_kde_2d( pc1_2[:,0], pc1_2[:,1], 
                                                        pc1_min, pc1_max, pc2_min, pc2_max, 
                                                        args.bandwidth,
                                                        xbins=args.xbins, ybins=args.ybins)

    pkl_data = {"pc1_grid":pc1_grid, "pc2_grid":pc2_grid, "density_grid":density_grid}

    out_pkl = args.out_prefix + "_" + name + ".pkl"
    pickle.dump( pkl_data, open(out_pkl, "w") )

    out_fig = args.out_prefix + "_" + name + ".pdf"
    plot_density(density_grid, pc1_min, pc1_max, pc2_min, pc2_max, "pc1", "pc2", out_fig)

