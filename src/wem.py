"""
The purpose of this script is to attempt to sample uniformly from the
set of all motifs having IC of c +/- epsilon.  We do this through the
weighted ensemble method.
"""
from random import shuffle,random,choice,randrange
from utils import *
from scipy.stats import norm
import mpmath as mp

def total_prob(binned_states):
    return sum(p for bs in binned_states for (state,p) in bs)

def downsample(states,M):
    p_total = sum(p for m,p in states)
    new_states = states[:]
    n = len(new_states)
    while n > M:
        i,j = randrange(n),randrange(n)
        if i == j:
            continue
        mi,pi = new_states[i]
        mj,pj = new_states[j]
        r = random.random()
        if r > pi/(pi+pj): #i gets killed
            new_states[j] = (mj,pi+pj)
            new_states.pop(i)
        else:              #j gets killed
            new_states[i] = (mi,pi+pj)
            new_states.pop(j)
        n = len(new_states)
    p_after = sum(p for m,p in new_states)
    return new_states

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
    
def upsample(states,M):
    deficiency = M - len(states)
    p_before = sum([p for state,p in states])
    ps = [p/p_before for m,p in states]
    added_states = [inverse_cdf_sample(states,ps)
                  for i in range(deficiency)]
    p_after = sum([p for state,p in added_states]) + p_before
    new_states = [(state,p*p_before/p_after)
                  for (state,p) in (states + added_states)]
    p_check = sum(map(second,new_states))
    return new_states

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


def weighted_ensemble(q, f, init_states, bins, M, tau, timesteps,final_bin_index=None,verbose=2):
    Q = lambda (state,p):(q(state),p)
    F = lambda (state,p):(f(state),p)
    n = float(len(init_states))
    n = mp.mpf(n)
    binned_states = [[(state,1/n)] for state in init_states]
    #print "total prob:",total_prob(binned_states)
    t = 0
    #history = []
    if final_bin_index:
        timesteps = 10**6
    while t < timesteps:
        next_tau = t + tau
        while t < next_tau:
            # Run the dynamics
            binned_states = mmap(Q,binned_states)
            t += 1
        new_binned_states = [[] for _ in pairs(bins)]
        for bs in binned_states:
            for state,p in bs:
                f_state = f(state)
                for i,(b0,b1) in enumerate(pairs(bins)):
                    if b0 <= f_state <b1:
                        index = i
                        break
                new_binned_states[i].append((state,p))
        binned_states = new_binned_states
        # trajectories are now resorted into correct bins
        for bs_index in range(len(binned_states)):
            binned_state = binned_states[bs_index]
            # If bin must be upsampled...
            if 0 < len(binned_state) < M:
                binned_states[bs_index] = upsample(binned_states[bs_index],M)
            #bin must be downsampled...
            elif len(binned_state) > M:
                binned_states[bs_index] = downsample(binned_states[bs_index],M)
        probs = map(lambda bs:sum(p for s,p in bs),binned_states)
        if verbose == 2:
            print t,probs
        elif verbose == 1:
            print t,sum(p > 0 for p in probs),len(bins)
        if final_bin_index and binned_states[final_bin_index]:
            return binned_states
    return binned_states

