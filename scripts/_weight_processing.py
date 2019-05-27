"""
functions to modify weights
"""
import copy

from _algdock import load_algdock_snapshots_for_each_of_six_yank_systems

def equalize_system_weights(original_weights):
    """
    make all the ref systems equally weighted
    :param original_weights: dict
    :return: new_weights, dict
    """
    new_weights = copy.deepcopy(original_weights)
    for system in new_weights["systems"].keys():
        new_weights["systems"][system] = 1.0
    return new_weights


def take_6_apo(original_weights):
    """
    only snapshots from apo of YANK states get nonzero weights
    :param original_weights: dict
    :return: new_weights
    """
    new_weights = equalize_system_weights(original_weights)
    ordered_snapshots = load_algdock_snapshots_for_each_of_six_yank_systems()

    for system in new_weights["systems"].keys():
        for i, snapshot in enumerate(ordered_snapshots[system]):
            if i < len(ordered_snapshots[system]) - 6:
                new_weights[system][snapshot] = 0.
            else:
                new_weights[system][snapshot] = 1.

    return new_weights


def take_12_near_apo(original_weights):
    """
    :param original_weights: dict
    :return: dict
    """
    new_weights = equalize_system_weights(original_weights)
    ordered_snapshots = load_algdock_snapshots_for_each_of_six_yank_systems()

    for system in new_weights["systems"].keys():
        for i, snapshot in enumerate(ordered_snapshots[system]):
            if i < len(ordered_snapshots[system]) - 12:
                new_weights[system][snapshot] = 0.
    return new_weights


def take_24_near_apo(original_weights):
    """
    :param original_weights: dict
    :return: dict
    """
    new_weights = equalize_system_weights(original_weights)
    ordered_snapshots = load_algdock_snapshots_for_each_of_six_yank_systems()

    for system in new_weights["systems"].keys():
        for i, snapshot in enumerate(ordered_snapshots[system]):
            if i < len(ordered_snapshots[system]) - 24:
                new_weights[system][snapshot] = 0.
    return new_weights


def take_6_holo(original_weights):
    """
    only snapshots from apo of YANK states get nonzero weights
    :param original_weights: dict
    :return: new_weights
    """
    new_weights = equalize_system_weights(original_weights)
    ordered_snapshots = load_algdock_snapshots_for_each_of_six_yank_systems()

    for system in new_weights["systems"].keys():
        for i, snapshot in enumerate(ordered_snapshots[system]):
            if i < 6:
                new_weights[system][snapshot] = 1.
            else:
                new_weights[system][snapshot] = 0.

    return new_weights


def take_12_near_holo(original_weights):
    """
    :param original_weights: dict
    :return: dict
    """
    new_weights = equalize_system_weights(original_weights)
    ordered_snapshots = load_algdock_snapshots_for_each_of_six_yank_systems()

    for system in new_weights["systems"].keys():
        for i, snapshot in enumerate(ordered_snapshots[system]):
            if i >= 12:
                new_weights[system][snapshot] = 0.

    return new_weights


def take_24_near_holo(original_weights):
    """
    :param original_weights: dict
    :return: dict
    """
    new_weights = equalize_system_weights(original_weights)
    ordered_snapshots = load_algdock_snapshots_for_each_of_six_yank_systems()

    for system in new_weights["systems"].keys():
        for i, snapshot in enumerate(ordered_snapshots[system]):
            if i >= 24:
                new_weights[system][snapshot] = 0.

    return new_weights
