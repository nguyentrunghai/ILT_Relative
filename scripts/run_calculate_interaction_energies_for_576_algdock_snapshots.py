
import argparse
import os
import sys

from _algdock import SIX_YANK_SYSTEMS, load_algdock_snapshots_for_each_of_six_yank_systems
from _process_yank_outputs import extract_complex_receptor_ligand_crds

sys.path.append("/home/tnguye46/opt/src/Bpmf_with_FFT")
from Md_OpenMM import openmm_energy

parser = argparse.ArgumentParser()
parser.add_argument('--yank_dir', type=str , default="/home/tnguye46/T4_Lysozyme/Yank")
parser.add_argument('--amber_dir', type=str , default="/home/tnguye46/T4_Lysozyme/Yank_snapshots/amber_for_6_ligands")

parser.add_argument('--repeats', type=int, default=3)
parser.add_argument('--nc_file', type=str, default='complex-implicit.nc')
args = parser.parse_args()

TEMPERATURE = 300.                                                                                                                               
KB = 8.3144621E-3/4.184  # kcal/mol/K
BETA = 1. / TEMPERATURE / KB

snapshots_for_each_system = load_algdock_snapshots_for_each_of_six_yank_systems()

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

    complex_pot_energies  = openmm_energy(complex_prmtop_file, complex_crds, "OpenMM_OBC2") * BETA       # in kT
    receptor_pot_energies = openmm_energy(receptor_prmtop_file, receptor_crds, "OpenMM_OBC2") * BETA
    ligand_pot_energies = openmm_energy(ligand_prmtop_file, ligand_crds, "OpenMM_OBC2") * BETA

    assert complex_pot_energies.shape[0] == receptor_pot_energies.shape[0] == ligand_pot_energies.shape[0], "energies terms do not have the same len"

    interaction_energies = complex_pot_energies - receptor_pot_energies - ligand_pot_energies
    snapshots = snapshots_for_each_system[yank_ligand]

    assert interaction_energies.shape[0] == len(snapshots), "snapshots and interaction_energies do not have the same len"

    out_file = yank_ligand + ".dat"
    with open(out_file, "w") as handle:
        handle.write( "# interaction energies of algock snapshots for %s\n" %yank_ligand )
        handle.write( "# in kT\n" )
        handle.write( "# algdock snapshot indices     inter eneegies  complex pot energies   receptor pot energies   ligand pot energies\n" )

        for snapshot, inter_e, complex_e, receptor_e, ligand_e in zip(snapshots, interaction_energies, complex_pot_energies, receptor_pot_energies, ligand_pot_energies ):
            handle.write("%10s %20.10e %20.10e %20.10e %20.10e\n" %(snapshot, inter_e, complex_e, receptor_e, ligand_e) )

