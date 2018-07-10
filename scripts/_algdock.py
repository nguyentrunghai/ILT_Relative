
import os
import numpy as np


SIX_YANK_SYSTEMS = ["lysozyme.active.A__ABJ", "benzene.A__AAA", "lysozyme.inactive.A__AAS", "1-methylpyrrole.A__AAA", "phenol.A__AAA", "p-xylene.A__AAA"]

SIX_YANK_SYSTEMS_NAMES = {}
SIX_YANK_SYSTEMS_NAMES["lysozyme.active.A__ABJ"]   = "n-hexylbenzene"
SIX_YANK_SYSTEMS_NAMES["benzene.A__AAA"]           = "benzene"
SIX_YANK_SYSTEMS_NAMES["lysozyme.inactive.A__AAS"] = "($\pm$)-camphor"
SIX_YANK_SYSTEMS_NAMES["1-methylpyrrole.A__AAA"]   = "methylpyrrole"
SIX_YANK_SYSTEMS_NAMES["phenol.A__AAA"]            = "phenol"
SIX_YANK_SYSTEMS_NAMES["p-xylene.A__AAA"]          = "p-xylene"

FOUR_ACTIVE_YANK_SYSTEMS = ["lysozyme.active.A__ABJ", "benzene.A__AAA", "1-methylpyrrole.A__AAA", "p-xylene.A__AAA"]

SNAPSHOT_LIST_DIR   = "/home/tnguye46/T4_Lysozyme/scripts"
SNAPSHOT_LIST_FILES = ["the_5th_288_snapshots_names.dat", "the_10th_288_snapshots_names.dat"]
SNAPSHOT_LIST_FILES = [os.path.join(SNAPSHOT_LIST_DIR, file) for file in SNAPSHOT_LIST_FILES]

def load_algdock_snapshots():
    """
    snapshot_list_files :   list of str
    return
        sorted list of Bing's notation snapshots
    """
    snapshots = []
    for filename in SNAPSHOT_LIST_FILES:
        snapshots.extend( open(filename, "r").read().split() )
    snapshots.sort(key = lambda i: int(i))
    return snapshots

def ligand_name_to_group_code(ligand_names):
    """
    ligand_names    :   list of str
    return 
        group_code[ligand_full_name]["group"]   -> group
        group_code[ligand_full_name]["code"]    -> code
    """
    assert type(ligand_names) == list, "ligand_names must be a list"
    group_code = {}
    for ligand in ligand_names:
        group_code[ligand] = {}
        group_code[ligand]["group"] = ligand[: -3]
        group_code[ligand]["code"] = ligand[-3 :]
    return group_code

def load_algdock_snapshots_for_each_of_six_yank_systems():
    """
    """
    snapshots = load_algdock_snapshots()
    total_number_of_snapshots = len(snapshots)
    number_of_systems = len(SIX_YANK_SYSTEMS)

    number_of_snapshots_per_system = total_number_of_snapshots / number_of_systems
    assert number_of_snapshots_per_system * number_of_systems == total_number_of_snapshots, "either total number of snapshots or systems is wrong"

    begins_ends = [ ( i*number_of_snapshots_per_system, (i+1)*number_of_snapshots_per_system) for i in range(number_of_systems)]
    
    snapshots_for_each_system = {}
    snapshots_for_each_system["all"] = snapshots
    for system, (begin, end) in zip(SIX_YANK_SYSTEMS, begins_ends):
        snapshots_for_each_system[system] = snapshots[begin : end]
    return snapshots_for_each_system

def Bing_snapshot_indices_for_algdock(nframes_per_state=10000):
    """
    """
    first_stride = 1000
    first_list = [i for i in range(nframes_per_state) if i%first_stride == 0]

    last_list = [first_list[5+10*i] for i in range(len(first_list)) if 5+10*i < len(first_list) ]
    last_list.extend( [ first_list[9+10*i] for i in range(len(first_list)) if 9+10*i < len(first_list) ] )
    return sorted(last_list)


def load_bpmfs(file_name, exclude_nan=False):
    """
    """
    bpmfs = {}
    with open(file_name, "r") as handle:
        for line in handle:
            snapshot, value = line.split()

            if value.lower() != "nan":
                bpmfs[snapshot] = np.float(value)
            else:
                if not exclude_nan:
                    bpmfs[snapshot] = np.nan
    return bpmfs


