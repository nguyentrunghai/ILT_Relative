
import argparse
import os
import pickle

import pandas as pd

from _kde_kullback_leibler_divergence import kullback_leibler_divergence, jensen_shannon_divergence, bin_area, annotated_heat_map

from _yank import YANK_LIGANDS

parser = argparse.ArgumentParser()

parser.add_argument( "--den_file_dir",  type=str, default="./")
parser.add_argument( "--out_kl_prefix",  type=str, default="kld_matrix")
parser.add_argument( "--out_js_prefix",  type=str, default="jsd_matrix")

args = parser.parse_args()

ref_ensembles = ["apo"] + ["1-methylpyrrole.A__AAA", "benzene.A__AAA", "p-xylene.A__AAA", 
                        "phenol.A__AAA", "lysozyme.inactive.A__AAS", "lysozyme.active.A__ABJ"] 

all_yank_24_holos = YANK_LIGANDS.keys()
target_18_holos = [ligand for ligand in all_yank_24_holos if ligand not in ref_ensembles]


kull_lei_d = {}
jen_shan_d = {}

kld_ligand_code = {}
jsd_ligand_code = {}

for ref_ensemble in ref_ensembles:
    
    ref_ensemble_name = "apo" if ref_ensemble == "apo" else YANK_LIGANDS[ref_ensemble]
    print "ref_ensemble_name: " + ref_ensemble_name

    kull_lei_d[ref_ensemble_name] = {}
    jen_shan_d[ref_ensemble_name] = {}

    kld_ligand_code[ref_ensemble] = {}
    jsd_ligand_code[ref_ensemble] = {}

    ref_den_file = os.path.join(args.den_file_dir, "density_" + ref_ensemble + ".pkl")
    print "loading " + ref_den_file + "\n"
    ref_data = pickle.load(open(ref_den_file, "r"))

    ref_density = ref_data["density_grid"]
    ref_bin_area = bin_area(ref_data["pc1_grid"], ref_data["pc2_grid"])

    for target_holo in target_18_holos:

        target_holo_name = YANK_LIGANDS[target_holo]
        print "target_holo_name: " + target_holo_name

        target_den_file = os.path.join(args.den_file_dir, "density_" + target_holo + ".pkl")
        print "loading " + target_den_file + "\n"
        target_data = pickle.load(open(target_den_file, "r"))

        target_density = target_data["density_grid"]
        target_bin_area = bin_area(target_data["pc1_grid"], target_data["pc2_grid"])

        assert target_bin_area == ref_bin_area, \
                "ref_bin_area %0.5f not the same as target_bin_area %0.5f" %(ref_bin_area, target_bin_area)

        kull_lei_d[ref_ensemble_name][target_holo_name] = kullback_leibler_divergence(target_density, ref_density, ref_bin_area)
        jen_shan_d[ref_ensemble_name][target_holo_name] = jensen_shannon_divergence(target_density, ref_density, ref_bin_area)

        kld_ligand_code[ref_ensemble][target_holo] = kull_lei_d[ref_ensemble_name][target_holo_name]
        jsd_ligand_code[ref_ensemble][target_holo] = jen_shan_d[ref_ensemble_name][target_holo_name]

pickle.dump( kull_lei_d, open(args.out_kl_prefix + ".pkl", "w") )
pickle.dump( jen_shan_d, open(args.out_js_prefix + ".pkl", "w") )

pickle.dump( kld_ligand_code, open(args.out_kl_prefix + "_ligand_code.pkl", "w") )
pickle.dump( jsd_ligand_code, open(args.out_js_prefix + "_ligand_code.pkl", "w") )


ref_ensemble_names = ["apo"] + [YANK_LIGANDS[ligand] for ligand in ref_ensembles if ligand != "apo"]
targ_holo_names = [YANK_LIGANDS[ligand] for ligand in target_18_holos]

kull_lei_dframe = pd.DataFrame(kull_lei_d, columns=ref_ensemble_names, index=targ_holo_names)

means_medians = pd.DataFrame( {"mean":kull_lei_dframe.mean(axis=0), "median":kull_lei_dframe.median(axis=0)} )
kull_lei_dframe = kull_lei_dframe.append(means_medians.T)

annotated_heat_map(kull_lei_dframe, args.out_kl_prefix + ".pdf" )



jen_shan_dframe = pd.DataFrame(jen_shan_d, columns=ref_ensemble_names, index=targ_holo_names)

means_medians = pd.DataFrame( {"mean":jen_shan_dframe.mean(axis=0), "median":jen_shan_dframe.median(axis=0)} )
jen_shan_dframe = jen_shan_dframe.append(means_medians.T)

annotated_heat_map(jen_shan_dframe, args.out_js_prefix + ".pdf" )

