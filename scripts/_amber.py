
import sys
sys.path.append("/home/tnguye46/opt/src/BPMFwFFT/bpmfwfft")
from IO import PrmtopLoad


def receptor_ligand_separation(complex_prmtop_file, receptor_prmtop_file, ligand_prmtop_file):
    """
    :param complex_prmtop_file: str
    :param receptor_prmtop_file: str
    :param ligand_prmtop_file: str
    :return: (receptor_indices, ligand_indices)
            receptor_indices: list of int, list of receptor atom index
            ligand_indices: list of int, list of ligand atom index
    """
    complex_prmtop = PrmtopLoad(complex_prmtop_file).get_all_parameters()
    receptor_prmtop = PrmtopLoad(receptor_prmtop_file).get_all_parameters()
    ligand_prmtop = PrmtopLoad(ligand_prmtop_file).get_all_parameters()

    assert len(ligand_prmtop["RESIDUE_LABEL"]) == 1, "There are none or more than one ligand molecules"
    assert len(ligand_prmtop["RESIDUE_LABEL"]) + len(receptor_prmtop["RESIDUE_LABEL"]) == len( complex_prmtop["RESIDUE_LABEL"] ), \
            "Number of residues of receptor plus ligand is not the same as complex"

    ligand_resname = ligand_prmtop["RESIDUE_LABEL"][0]

    n_complex_atoms = complex_prmtop["POINTERS"]["NATOM"]
    n_recptor_atoms = receptor_prmtop["POINTERS"]["NATOM"]
    n_ligand_atoms = ligand_prmtop["POINTERS"]["NATOM"]

    if complex_prmtop["RESIDUE_LABEL"][0] == ligand_resname:
        ligand_indices = range(0, n_ligand_atoms)
        receptor_indices = range(n_ligand_atoms, n_complex_atoms)
        
    elif complex_prmtop["RESIDUE_LABEL"][-1] == ligand_resname:
        receptor_indices = range(0, n_recptor_atoms)
        ligand_indices = range(n_recptor_atoms, n_complex_atoms)

    else:
        raise Exception("Ligand is not either on top or bottom of complex")

    return receptor_indices, ligand_indices


