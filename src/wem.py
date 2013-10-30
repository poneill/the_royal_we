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

def random_motif_with_dirty_bits(length,num_sites):
    motif = random_motif(length,num_sites)
    dirty_bits = [False] * length
    ics = map(lambda col:2-dna_entropy(col),
              transpose(motif))
    return motif,ics

def mutate_motif_with_dirty_bits((motif,ics)):
    length = len(motif[0])
    num_sites = len(motif)
    i = random.randrange(num_sites)
    j = random.randrange(length)
    site = motif[i]
    b = site[j]
    new_b = random.choice([c for c in "ACGT" if not c == b])
    new_site = subst(site,new_b,j)
    new_motif = [site if k != i else new_site
                 for k,site in enumerate(motif)]
    jth_column = transpose(new_motif)[j]
    jth_ic = 2-dna_entropy(jth_column)
    new_ics = subst(ics,jth_ic,j)
    return new_motif,new_ics

def motif_ic_with_dirty_bits((motif,ics)):
    return sum(ics)

def test_dirty_bits():
    for i in range(1000):
        m = random_motif_with_dirty_bits(10,16)
        m_ic1 = motif_ic_with_dirty_bits(m)
        m_ic2 = motif_ic(m[0])
        m1 = mutate_motif_with_dirty_bits(m)
        m1_ic1 = motif_ic_with_dirty_bits(m1)
        m1_ic2 = motif_ic(m1[0])
        print m_ic1 - m_ic2
        print m1_ic1 - m1_ic2
        
def sample_motifs(length,num_sites,ic,epsilon):
    bins = [-10] + myrange(0,ic,epsilon) + [ic,ic + epsilon]+ [ic + 10]
    tau = 1
    timesteps = 50
    M = 100
    init_states = [random_motif(10,16) for i in range(M)]
    results = weighted_ensemble(mutate_motif, motif_ic, init_states,
                                bins, M, tau, timesteps)
    motifs,ps = transpose(concat(results[-3:-1]))
    return inverse_cdf_sample(motifs,normalize(ps))

def sample_motifs_with_dirty_bits(length,num_sites,ic,epsilon,inv_cdf=True,timesteps = 100,wait_indefinitely=False,init_states=None):
    print "timesteps:",timesteps
    bins = [-10] + myrange(0,ic,epsilon) + [ic,ic + epsilon]+ [ic + 10]
    final_bin_index = None if not wait_indefinitely else -2
    tau = 1
    M = 100
    if init_states == None:
        init_states = [random_motif_with_dirty_bits(length,num_sites)
                       for i in range(M)]
    results = weighted_ensemble(mutate_motif_with_dirty_bits,
                                motif_ic_with_dirty_bits,
                                init_states,
                                bins, M, tau, timesteps,
                                final_bin_index=final_bin_index,verbose=1)
    try:
        motifs,ps = transpose(concat(results[-3:-1]))
    except ValueError:
        motifs,ps = transpose(concat(results))
        return sample_motifs_with_dirty_bits(length,num_sites,ic,epsilon,inv_cdf,timesteps, wait_indefinitely,init_states=motifs)
    if inv_cdf:
        motif, ics = inverse_cdf_sample(motifs,normalize(ps))
        return motif
    else:
        stripped_motifs = [motif for motif,ics in motifs]
        return stripped_motifs

def match_motif(motif,epsilon=0.1,inv_cdf=True,wait_indefinitely=False):
    length = len(motif[0])
    num_sites = len(motif)
    ic = motif_ic(motif)
    return sample_motifs_with_dirty_bits(length,num_sites,ic,epsilon,inv_cdf=inv_cdf,
                                         wait_indefinitely=wait_indefinitely)
    
def wem_control_motif(motif,epsilon=0.1):
    num_sites = len(motif)
    length = len(motif[0])
    ic = motif_ic(motif)
    return sample_motifs(length,num_sites,ic,epsilon)

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

def ic_we():
    bins = range(-1,3) + myrange(3,12,0.1) + [50]
    hist = weighted_ensemble(mutate_motif,
                             motif_ic,
                             [random_motif(10,16) for i in range(100)],
                             bins=bins,
                             M=1000,
                             tau=1,
                             timesteps=500)
    return hist
    
print "loaded"
