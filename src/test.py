"""
This file contains code for testing wem.py
"""

def test_we():
    """Check the validity of the implementation by sampling the random
    walk on the 1d integer lattice"""
    min_bin = 0
    max_bin = 100
    def q(x):
        if x == min_bin:
            return min_bin + 1
        elif x == max_bin:
            return max_bin - 1
        else:
            return x- 1 if random.random() < .5 else x + 1
    f = lambda x:x
    init_states = [0]
    bins = range(0,max_bin,max_bin/10) + [100]
    M = 1000
    tau = 1000
    timesteps = tau * 100
    results = weighted_ensemble(q, f, init_states, bins, M, tau, timesteps)
    return results
