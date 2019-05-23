
from __future__ import print_function

import os
import numpy as np
import netCDF4 as nc

from _algdock import Bing_snapshot_indices_for_algdock
from _algdock import SIX_YANK_SYSTEMS
from _amber import receptor_ligand_separation


def read_reduced_energy(yank_nc_files):
    """
    :param yank_nc_files: list of str, YANK *.nc output files
    :return: (u_kln, N_k)
              u_kln: 3d numpy array,  reduced potential energy
              N_k: 1d numpy array, number of samples at state k
    """
    assert len(yank_nc_files) > 1, "no nc file to read"
    nc_handles = [ nc.Dataset(filename, "r") for filename in yank_nc_files ]

    nsamples_per_state = np.array( [ handle.variables["positions"].shape[0] for handle in nc_handles ], dtype=int )
    total_nsamples = nsamples_per_state.sum()

    nstates = np.array( [ handle.variables["positions"].shape[1] for handle in nc_handles ], dtype=int  )
    if not np.all( nstates == nstates.mean() ):
        raise Exception("number of states K are not the same in different nc files")

    K = nstates[0]
    N_k = np.ones([K], dtype=int)*total_nsamples
    u_kln = np.zeros( [ K, K, N_k.max() ], dtype=float)
    #
    count = -1
    for handle in nc_handles:
        energies = handle.variables["energies"]
        nframes  = handle.variables["positions"].shape[0]

        for frame in range(nframes):
            count += 1
            state_indices = handle.variables["states"][frame, :]
            u_kln[state_indices, :, count] = energies[frame, :, :] 
    #
    return u_kln, N_k


def read_reduced_energies_for_algdock_snapshots(yank_nc_files):
    """
    :param yank_nc_files: list of str, YANK *.nc output files
    :return: (u_kln, N_k)
             u_kln: 3d numpy array,  reduced potential energy
             N_k: 1d numpy array, number of samples at state k
    """
    assert len(yank_nc_files) > 1, "no nc file to read"
    nc_handles = [ nc.Dataset(filename, "r") for filename in yank_nc_files ]

    seleted_frames_per_state = Bing_snapshot_indices_for_algdock(nframes_per_state=10000)

    nstates = np.array( [ handle.variables["positions"].shape[1] for handle in nc_handles ], dtype=int  )
    if not np.all( nstates == nstates.mean() ):
        raise Exception("number of states K are not the same in different nc files")

    nstates = nstates[0]
    nsamples_per_state = len(seleted_frames_per_state) * len(nc_handles)
    print("nsamples_per_state", nsamples_per_state)

    K = nstates
    N_k = np.ones([K], dtype=int) * nsamples_per_state
    u_kln = np.zeros( [ K, K, N_k.max() ], dtype=float)

    count = -1
    for handle in nc_handles:
        energies = handle.variables["energies"]
        for frame in seleted_frames_per_state:
            count += 1
            state_indices = handle.variables["states"][frame, :]
            u_kln[state_indices, :, count] = energies[frame, :, :]

    return u_kln, N_k


def augmented_energy_matrix_for_holo_marginal_weights( u_kln, N_k,
                                                        bpmfs, 
                                                        receptor_pot_energies):
    """
    Add a state into the reduced potential energy matrix, u_kln[K, K, N] ->  aug_u_kln[K, K+1, N]
    This new state is the marginal holo, with the ligand degrees of freedom being integrated out
    :param u_kln: 3d numpy array,  reduced potential energy
    :param N_k: 1d numpy array, number of samples at state k
    :param bpmfs: numpy array
    :param receptor_pot_energies:  numpy array
    :return: (aug_u_kln, aug_N_k)
            aug_u_kln: 3d numpy array,  reduced potential energy
            aug_N_k: number of samples at state k
    """
    assert receptor_pot_energies.shape == bpmfs.shape, "receptor_pot_energies and bpmfs must have the same len"
    K, _, N = u_kln.shape
    aug_N_k = np.zeros([K+1], dtype=int)
    aug_N_k[:K] = N_k[:]
    aug_u_kln = np.zeros( [ K, K+1, N ], dtype=float)

    aug_u_kln[:, 0:K, :] = u_kln[:, :, :]

    bpmfs[ np.isnan(bpmfs) ] = np.inf
    aug_u_kln[:, K, :] = (receptor_pot_energies + bpmfs).reshape(K, N)

    return aug_u_kln, aug_N_k 


def extract_reduced_interaction_energies_for_holo(u_kln):
    """
    :param u_kln:  3d numpy array,  reduced potential energy
    :return: interaction_energies, list of float
    """
    seleted_frames = Bing_snapshot_indices_for_algdock(nframes_per_state=10000*3)

    nstates = u_kln.shape[0]

    interaction_energies = []
    for state in range(nstates):
        for frame in seleted_frames:
            interaction_energies.append( u_kln[state, 0, frame] - u_kln[state, -1, frame] )
    return interaction_energies


def extract_complex_receptor_ligand_crds(yank_nc_files, complex_prmtop_file,
                                         receptor_prmtop_file, ligand_prmtop_file):
    """
    extract coordinated for Bing's snapshots

    :param yank_nc_files: list of YANK *.nc output files
    :param complex_prmtop_file: str
    :param receptor_prmtop_file: str
    :param ligand_prmtop_file: str
    :return: (complex_crds, receptor_crds, ligand_crds)
    """
    assert len(yank_nc_files) > 1, "no nc file to read"
    receptor_indices, ligand_indices = receptor_ligand_separation(complex_prmtop_file, receptor_prmtop_file, ligand_prmtop_file)
    selected_frames = Bing_snapshot_indices_for_algdock(nframes_per_state=10000)

    nc_handles = [ nc.Dataset(filename, "r") for filename in yank_nc_files ]
    nsamples_per_state = np.array( [ handle.variables["positions"].shape[0] for handle in nc_handles ], dtype=int )
    #total_nsamples = nsamples_per_state.sum()

    nstates = np.array( [ handle.variables["positions"].shape[1] for handle in nc_handles ], dtype=int  )
    if not np.all( nstates == nstates.mean() ):
        raise Exception("number of states K are not the same in different nc files")

    nstates = nstates[0]
    count = -1
    complex_crds = []
    receptor_crds = []
    ligand_crds = []
    for state in range(nstates):

        for handle in nc_handles:

            for frame in selected_frames:

                count += 1
                state_indices = handle.variables["states"][frame, :]
                replica_index = list(state_indices).index(state)

                complex_crd  = handle.variables['positions'][frame, replica_index, :, : ]*10.0
                receptor_crd = handle.variables['positions'][frame, replica_index, receptor_indices, : ]*10.0
                ligand_crd   = handle.variables['positions'][frame, replica_index, ligand_indices, : ]*10.0

                complex_crds.append(complex_crd)
                receptor_crds.append(receptor_crd)
                ligand_crds.append(ligand_crd)
    
    return np.array(complex_crds), np.array(receptor_crds), np.array(ligand_crds)


def load_interaction_energies(snaphot_col=0, value_col=1,
                                path="OpenMM_OBC2_interaction_energies_for_576_algdock_snapshots"):
    """
    :param snaphot_col: int
    :param value_col: int
    :param path: str
    :return: interaction_energies, dict
            interaction_energies[ligand][snapshot] -> float
    """
    interaction_energies = {}
    for ligand in SIX_YANK_SYSTEMS:
        filename = os.path.join(path, ligand+".dat")
        print("loaing "+filename)
        
        interaction_energies[ligand] = {}
        with open(filename, "r") as handle:
            for line in handle:
                if "#" not in line:
                    entries = line.split()
                    snapshot = entries[snaphot_col]
                    value = entries[value_col]
                    interaction_energies[ligand][snapshot] = np.float(value)
    return interaction_energies


def load_receptor_pot_energies(file_name, value_col=3):
    """
    :param file_name: str
    :param value_col: int
    :return: energies, numpy array
    """
    energies = []
    with open(file_name, "r") as handle:
        for line in handle:
            if "#" not in line:
                entries = line.split()
                value    = entries[value_col]
                energies.append(np.float(value))
    return np.array(energies)

