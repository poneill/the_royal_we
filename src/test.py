"""
This file contains code for testing wem.py
"""
import random
from wem import *

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

def test_downsample():
    before_num = 2000
    after_num = 1000
    xs = [random.normalvariate(0,1) for i in range(before_num)]
    ps = [1/float(before_num) for i in range(before_num)]
    states = zip(xs,ps)
    mu0 = sum([x*p for x,p in states])
    sigma0 = sum([(x**2)*p for x,p in states])
    new_states = downsample(states,after_num)
    mu1  = sum([x*p for x,p in new_states])
    sigma1 = sum([(x**2)*p for x,p in new_states])
    print mu0,sigma0
    print mu1,sigma1

def test_upsample():
    before_num = 500
    after_num = 1000
    xs = [random.normalvariate(0,1) for i in range(before_num)]
    ps = [1/float(before_num) for i in range(before_num)]
    states = zip(xs,ps)
    mu0 = sum([x*p for x,p in states])
    sigma0 = sum([(x**2)*p for x,p in states])
    new_states = upsample(states,after_num)
    mu1  = sum([x*p for x,p in new_states])
    sigma1 = sum([(x**2)*p for x,p in new_states])
    print mu0,sigma0
    print mu1,sigma1
