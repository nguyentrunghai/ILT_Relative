
from __future__ import print_function

import numpy as np
import os
import pickle

from _algdock import Bing_snapshot_indices_for_algdock

WEIGHTS_DIR = "/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/Mbar_Weights_holo_OBC2"

SNAPSHOT_LIST_DIR = "/home/tnguye46/T4_Lysozyme/Relative_Binding_FE/snapshot_list"
SNAPSHOT_LIST_FILES = ["the_5th_288_snapshots_names.dat", "the_10th_288_snapshots_names.dat"]
SNAPSHOT_LIST_FILES = [os.path.join(SNAPSHOT_LIST_DIR, file) for file in SNAPSHOT_LIST_FILES]

#$$$$$$$$$$$$$$$
ASSIGNMENT_FILE = "/home/tnguye46/T4_Lysozyme/Yank_snapshots/assign_snapshots_to_576_classes/assign.pkl"
ASSIGNMENT = pickle.load( open(ASSIGNMENT_FILE, "r") )

DOCK_STRIDE = 5000
FIRST_INDEX = 4999

NSTATES = 16
NR_YANK_SNAPSHOTS_PER_STATE = 10000
NR_YANK_SNAPSHOTS_PER_SYSTEM = 3 * 16 * NR_YANK_SNAPSHOTS_PER_STATE
from _algdock import SIX_YANK_SYSTEMS as YANK_SYSTEMS 


def load_snapshot_list():
    snapshots = []
    for file in SNAPSHOT_LIST_FILES:
        snapshots.extend( open(file, "r").read().split() )
    snapshots.sort(key = lambda i: int(i))
    return snapshots


def load_mbar_weights():
    snapshots = load_snapshot_list()
    nr_systems = len(YANK_SYSTEMS)

    nr_of_states_per_system = NR_YANK_SNAPSHOTS_PER_SYSTEM // NR_YANK_SNAPSHOTS_PER_STATE
    print("number of states per system %d"%nr_of_states_per_system)
    
    nr_selected_snapshots_per_system = NR_YANK_SNAPSHOTS_PER_SYSTEM // DOCK_STRIDE
    print('number of selected snapshots %d' %(nr_selected_snapshots_per_system * nr_systems))

    #selected_snapshots_per_system = range(FIRST_INDEX, NR_YANK_SNAPSHOTS_PER_SYSTEM, DOCK_STRIDE)
    selected_snapshots_per_system = Bing_snapshot_indices_for_algdock(nframes_per_state=NR_YANK_SNAPSHOTS_PER_SYSTEM)


    if len(selected_snapshots_per_system) != nr_selected_snapshots_per_system:
        raise RuntimeError("len(selected_snapshots_per_system) != nr_selected_snapshots_per_system")

    if len(snapshots) != (nr_selected_snapshots_per_system * nr_systems):
        raise RuntimeError("nr of snapshot from the list file is %d but number estimated from stride is %d "%(
            len(snapshots), nr_selected_snapshots_per_system * nr_systems))

    print 'parsing MBAR weights ...'
    #
    bl_weights = {}
    st_weights = {}
    single_weights = {}

    system_wide_weights = {}

    index_range_in_each_state = [(i*NR_YANK_SNAPSHOTS_PER_STATE, (i+1)*NR_YANK_SNAPSHOTS_PER_STATE) for i in range(nr_of_states_per_system) ]
    repeated_index_range_in_each_state = []
    ## repeated twice, not a bug
    for ind in index_range_in_each_state:
        repeated_index_range_in_each_state.append(ind)
        repeated_index_range_in_each_state.append(ind)      # repeated twice, not a bug
    print repeated_index_range_in_each_state

    #$$$$$$$$$$$$$$$$$$$$$$$$$$
    all_weights = {}

    for system in YANK_SYSTEMS:
        weights = np.loadtxt( os.path.join(WEIGHTS_DIR, system + ".weig") )[:, -1]     
        system_wide_weights[system] = weights[-len(weights)/NSTATES : ].sum()
        
        bl_weights[system] = np.array( [ weights[i*DOCK_STRIDE : (i+1)*DOCK_STRIDE].sum() for i in range(len(selected_snapshots_per_system)) ] )
        bl_weights[system] /= bl_weights[system].sum()

        st_weights[system]  = np.array( [ weights[i : j].sum() for i, j in repeated_index_range_in_each_state ] )
        st_weights[system]  /= st_weights[system].sum()

        single_weights[system] = np.array( [ weights[s] for s in selected_snapshots_per_system ] )
        single_weights[system] /= single_weights[system].sum()

        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        all_weights[system] = weights

    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    system_wide_weights_normalized = {}
    for system in YANK_SYSTEMS:
        system_wide_weights_normalized[system] = system_wide_weights[system] / np.sum(system_wide_weights.values())

    all_weights_equal_sys   = np.concatenate( [ all_weights[system] for system in YANK_SYSTEMS ] )
    all_weights_ub_weighted = np.concatenate( [ all_weights[system]*system_wide_weights_normalized[system] for system in YANK_SYSTEMS ] )

    # returned dicts
    block_weights = {"systems" : system_wide_weights} # mbar_weights[system][snapshot]
    state_weights = {"systems" : system_wide_weights}
    single_snap_weights = {"systems" : system_wide_weights}

    #$$$$$$$$$$$$$$$$$$$$$$$$$$
    stru_group_weights_equal_sys   = {"systems" : {system:1./len(YANK_SYSTEMS) for system in YANK_SYSTEMS} }
    stru_group_weights_ub_weighted = {"systems" : {system:1./len(YANK_SYSTEMS) for system in YANK_SYSTEMS} }

    block_w = {}
    count = -1
    for system in YANK_SYSTEMS:
        for w in bl_weights[system]:
            count += 1
            block_w[ snapshots[count] ] = w

    state_w = {}
    count = -1
    for system in YANK_SYSTEMS:
        for w in st_weights[system]:
            count += 1
            state_w[ snapshots[count] ] = w

    single_snap_w = {}
    count = -1
    for system in YANK_SYSTEMS:
        for w in single_weights[system]:
            count += 1
            single_snap_w[ snapshots[count] ] = w

    #$$$$$
    stru_group_w_equal_sys   = {}
    stru_group_w_ub_weighted = {}
    for snapshot in ASSIGNMENT.keys():
        stru_group_w_equal_sys[snapshot]   = all_weights_equal_sys[ np.array(ASSIGNMENT[snapshot], dtype=int) ].sum()
        stru_group_w_ub_weighted[snapshot] = all_weights_ub_weighted[ np.array(ASSIGNMENT[snapshot], dtype=int) ].sum() 

    #---------------------------------------------------
    snapshot_for_each_system = {}
    count = -1
    for system in YANK_SYSTEMS:
        snapshot_for_each_system[system] = []
        for i in range(nr_selected_snapshots_per_system):
            count += 1
            snapshot_for_each_system[system].append( snapshots[count] )
    #
    for system in YANK_SYSTEMS:
        block_weights[system] = {}
        for snapshot in snapshots:
            if snapshot in snapshot_for_each_system[system]:
                block_weights[system][snapshot] = block_w[snapshot]

    for system in YANK_SYSTEMS:
        state_weights[system] = {}
        for snapshot in snapshots:
            if snapshot in snapshot_for_each_system[system]:
                state_weights[system][snapshot] = state_w[snapshot]

    for system in YANK_SYSTEMS:
        single_snap_weights[system] = {}
        for snapshot in snapshots:
            if snapshot in snapshot_for_each_system[system]:
                single_snap_weights[system][snapshot] = single_snap_w[snapshot]

    #$$$$$$$$$$$$$$$$$$$$$$
    for system in YANK_SYSTEMS:
        stru_group_weights_equal_sys[system] = {}
        stru_group_weights_ub_weighted[system] = {}

        for snapshot in snapshots:
            if snapshot in snapshot_for_each_system[system]:

                stru_group_weights_equal_sys[system][snapshot]   = stru_group_w_equal_sys[snapshot]
                stru_group_weights_ub_weighted[system][snapshot] = stru_group_w_ub_weighted[snapshot]

    #
    return block_weights, state_weights, single_snap_weights, stru_group_weights_equal_sys, stru_group_weights_ub_weighted




