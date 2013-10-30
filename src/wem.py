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

def weighted_ensemble(q, f, init_states, bins, M, tau, timesteps,
                      final_bin_index=None,verbose=2):
    """
    Sample the distribution associated with the Markov process given
    by transition kernel (or matrix) q by applying the Weighted
    Ensemble Method along the coordinate specified by bins.

    Input parameters:

    q - a Markov kernel (matrix)
    
    f - a projection function mapping states onto the reaction coordinate (rc)

    init_states - a set of initial states for seeding the ensemble

    bins - a list specifying the left endpoints of the bins along the rc.

    M - number of states per occupied bin.

    timesteps - number of iterations of q.

    tau - number of iterations of q per resorting of bins.

    final_bin_index - if set to non-negative integer j, halt
    trajectories and return as soon as the jth bin becomes occupied.

    verbose - set verbosity (in interval [0,2]).
    """
    # The basic datatype is a state associated with a probability.  We
    # need q and f to act on those types, so define:
    Q = lambda (state,p):(q(state),p)
    F = lambda (state,p):f(state)
    n = float(len(init_states))
    # We will also be dealing with very low probabilities, so use
    # multiple precision here.  This turns out to have a negligible
    # effect on performance.
    n = mp.mpf(n)
    # assign init_states uniform probabilities initially
    binned_states = [[(state,1/n)] for state in init_states]
    t = 0
    while t < timesteps:
        next_tau = t + tau
        while t < next_tau:
            # Run the dynamics, iterating every state with Q.
            binned_states = mmap(Q,binned_states)
            t += 1

        # Now resort the states into the correct bins.  We do this in
        # an ugly imperative fashion for performance reasons, because
        # Python likes to punish you for writing pretty code.
        new_binned_states = [[] for _ in pairs(bins)]
        for bs in binned_states:
            for state_p in bs:
                # (1) project state onto reaction coord
                f_state = F(state_p)
                # (2) find the right bin for it
                for i,(b0,b1) in enumerate(pairs(bins)):
                    if b0 <= f_state <b1:
                        index = i
                        break
                # (3) and put the state in that bin
                new_binned_states[i].append(state_p)
        binned_states = new_binned_states
        # trajectories are now resorted into correct bins

        # But now some bins may contain more than M states, some
        # less...

        # We also do things here in a slightly ugly fashion to deal
        # with the fact that we are mutating a list while iterating
        # over it, which is despised in Python.
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

