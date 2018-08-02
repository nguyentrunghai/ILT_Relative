
import argparse
import os
import sys

from _algdock import SIX_YANK_SYSTEMS, load_algdock_snapshots_for_each_of_six_yank_systems
from _process_yank_outputs import extract_complex_receptor_ligand_crds

sys.path.append("/home/tnguye46/opt/src/Bpmf_with_FFT")
from Md_OpenMM import openmm_energy
import IO

parser = argparse.ArgumentParser()
parser.add_argument('--yank_dir', type=str , default="/home/tnguye46/T4_Lysozyme/Yank")
parser.add_argument('--amber_dir', type=str , default="/home/tnguye46/T4_Lysozyme/Yank_snapshots/amber_for_6_ligands")

parser.add_argument('--repeats', type=int, default=3)
parser.add_argument('--nc_file', type=str, default='complex-implicit.nc')
args = parser.parse_args()

snapshots_for_each_system = load_algdock_snapshots_for_each_of_six_yank_systems()
print snapshots_for_each_system.keys()

for yank_ligand in SIX_YANK_SYSTEMS:
    yank_nc_files = [ os.path.join(args.yank_dir, yank_ligand, "Repeat%d"%repeat, "output", args.nc_file) 
                            for repeat in range(args.repeats) ]

    complex_prmtop_file  = os.path.join(args.amber_dir, yank_ligand, "complex.prmtop")
    receptor_prmtop_file = os.path.join(args.amber_dir, yank_ligand, "receptor.prmtop")
    ligand_prmtop_file   = os.path.join(args.amber_dir, yank_ligand, "ligand.prmtop")

    print "\n".join(yank_nc_files)
    print complex_prmtop_file
    print receptor_prmtop_file
    print ligand_prmtop_file
    print "\n"

    complex_crds, receptor_crds, ligand_crds = extract_complex_receptor_ligand_crds(yank_nc_files, complex_prmtop_file, receptor_prmtop_file, ligand_prmtop_file)

    snapshots = snapshots_for_each_system[yank_ligand]

    assert  receptor_crds.shape[0] == len(snapshots), "snapshots and  receptor_crds do not have the same len"

    for i, snapshot in enumerate(snapshots):
        complex_out = "complex_" +  snapshot + ".pdb"
        IO.write_pdb(complex_prmtop_file, complex_crds[i], complex_out, "w")

        receptor_out = "receptor_" +  snapshot + ".pdb"
        IO.write_pdb(receptor_prmtop_file, receptor_crds[i], receptor_out, "w")

        ligand_out = "ligand_" +  snapshot + ".pdb"
        IO.write_pdb(ligand_prmtop_file, ligand_crds[i], ligand_out, "w")

print "DONE"
