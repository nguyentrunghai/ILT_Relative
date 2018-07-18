
import numpy as np
import pymbar

def run_mbar(u_kln, N_k):
    """
    u_kln   :   ndarray, reduced energy matrix
    N_k     :   ndarray, number of samples at state k
    """
    K = len(N_k)
    f_k_BAR = np.zeros(K)
    for k in range(K-2):
        w_F = u_kln[ k, k+1, :N_k[k] ] - u_kln[ k, k, :N_k[k] ]
        w_R = u_kln[ k+1, k, :N_k[k+1] ] - u_kln[ k+1, k+1, :N_k[k+1] ]
        f_k_BAR[k+1] = pymbar.BAR(w_F, w_R, relative_tolerance=0.000001, \
                verbose=False, compute_uncertainty=False)
    f_k_BAR = np.cumsum(f_k_BAR)
    mbar = pymbar.MBAR(u_kln, N_k, verbose = True, initial_f_k = f_k_BAR)
    return mbar

def mbar_weights_in_state(u_kln, N_k, state_index):
    """
    u_kln   :   ndarray, reduced energy matrix
    N_k     :   ndarray, number of samples at state k
    state_index   :   int, usually for YANK 0 is fully interacting, K-1 is fully noninteracting
    """
    if state_index >= 0:
        assert 0 <= state_index < u_kln.shape[0], "state index out of range"
    else:
        assert 1 <= -state_index <= u_kln.shape[0], "state index out of range"
    mbar = run_mbar(u_kln, N_k)
    weights = mbar.getWeights()
    weights = weights[:, state_index]
    return weights

def mbar_weights_in_holo(u_kln, N_k):
    return mbar_weights_in_state(u_kln, N_k, 0)
