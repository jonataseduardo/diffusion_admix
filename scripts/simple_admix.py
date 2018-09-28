import moments
import matplotlib.pyplot as pl
import numpy as np
import numpy
import OutOfAfrica as ooa
import pandas
from collections import OrderedDict
from itertools import combinations



def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def onepop_moments(fs, k):
    ns = fs.sample_sizes[0]
    # sample frequencies p 
    p = numpy.arange(0, ns + 1, dtype=float) / ns
    return (fs *(p ** k)).sum()

def cross_moments(fs, k):
    r = fs.Npop
    ns = fs.sample_sizes

    # counts_per_pop is an r+1 dimensional array, where the last axis simply
    # records the indices of the entry. 
    # For example, counts_per_pop[4,19,8] = [4,19,8]
    counts_per_pop = numpy.indices(fs.shape)
    counts_per_pop = numpy.transpose(counts_per_pop, axes=range(1, r + 1) + [0])
    
    ktilde = np.tile(k, np.append(ns + 1, 1))

    x = 1. * counts_per_pop / ns 

    x_k = numpy.power(x, ktilde)
    
    # Included tuple() to avoid the following warning msg
    # FutureWarning: Using a non-tuple sequence for multidimensional 
    # indexing is deprecated
    this_slice = tuple([slice(None)] * r + [numpy.newaxis])

    fstilde = numpy.tile(fs.data[this_slice], [1] * r  + [r])

    mu_k = x_k * fstilde

    return mu_k.sum()

def OutOfAfrica_stepwise((nuAf, nuB, nuEu0, nuEu, nuAs0, nuAs, 
                           mAfB, mAfEu, mAfAs, mEuAs, TAf, TB, TEuAs), 
                          (n1,n2,n3),
                          n_steps = None):
    """
    A three-population model used to model out-of-Africa demography.
    """

    Ttotal = TAf + TB + TEuAs

    niterAf = ((TAf / Ttotal) * n_steps).astype(np.integer)
    niterB = ((TB / Ttotal) * n_steps).astype(np.integer)
    niterEuAs = ((TEuAs / Ttotal) * n_steps).astype(np.integer)

    dtAf = TAf / niterAf
    dtB = TB / niterB
    dtEuAs = TEuAs / niterEuAs

    fs_history = OrderedDict()

    #first step: a single population
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2+n3)
    fs = moments.Spectrum(sts)
    
    dt = dtAf 
    time = 0
    step = 0 
    split_1 = split_2 = True
    while(time <= 1.001 * Ttotal ):
        time += dt

        if(time <= 1.001 * TAf):
            print('before split' , step, time )
            #integrate for time TAf (with constant population)
            fs.integrate([nuAf], dt, dt)
            fs_history[step] = [time, fs]
            step += 1

        if(time > TAf and time <= 1.001 * (TAf + TB)):
            if(split_1):
                print('first split' , step, time )
                split_1 = False
                dt = dtB
                fs = moments.Manips.split_1D_to_2D(fs, n1, n2+n3)
                fs_history[step] = [time, fs]
                step += 1
                mig1=numpy.array([[0, mAfB],[mAfB, 0]])

            print('after first split' , time )
            fs.integrate([nuAf, nuB], dt, 0.5 * dt, m=mig1)
            fs_history[step] = [time, fs]
            step += 1

        if(time > TAf + TB and time <= 1.001 * (TAf + TB + TEuAs)):
            if(split_2):
                print('sencond split' , step, time )
                split_2 = False
                dt = dtEuAs
                fs = moments.Manips.split_2D_to_3D_2(fs, n2, n3)
                fs_history[step] = [time, fs]
                step += 1
                # migration rates matrix
                mig2=numpy.array([[0, mAfEu, mAfAs],
                                  [mAfEu, 0, mEuAs],
                                  [mAfAs, mEuAs, 0]])
                nuEu0t, nuAs0t = nuEu0, nuAs0

            print('after sencond split' , step, time )

            #define functions for population sizes
            nuEu_func = lambda t: nuEu0t*(nuEu/nuEu0t)**(t/TEuAs)
            nuAs_func = lambda t: nuAs0t*(nuAs/nuAs0t)**(t/TEuAs)
            nu2 = lambda t: [nuAf, nuEu_func(t), nuAs_func(t)]

            # migration rates matrix
            mig2=numpy.array([[0, mAfEu, mAfAs],[mAfEu, 0, mEuAs],[mAfAs, mEuAs, 0]])
            fs.integrate(nu2, dt, 0.5 * dt, m=mig2)

            fs_history[step] = [time, fs]
            nuEu0t = nuEu0t*(nuEu/nuEu0t)**(dt/TEuAs)
            nuAs0t = nuAs0t*(nuAs/nuAs0t)**(dt/TEuAs)
            step += 1
                                
    return fs_history

def OOA_stats(fs_history, s, h):

    stats_pop = dict()
    stats_pop["step"] = []
    stats_pop["time"] = []
    stats_pop["pop_id"] = []
    stats_pop["pi"] = []
    stats_pop["load"] = []
    stats_pop["efficacy"] = []

    for step in fs_history.keys():
        np = fs_history[step][1].Npop
        rnp = range(np)
        lpop = [(list(set(rnp) ^ set(i)), i) for i in combinations(rnp, np - 1)]
        for i in lpop:
            stats_pop["step"] += [step]
            stats_pop["time"] += [fs_history[step][0]]
            stats_pop["pop_id"] += [i[0][0]]
            try:
                fs_i = fs_history[step][1].marginalize(i[1])
            except:
                fs_i = fs_history[step][1]

            stats_pop["pi"] += [fs_i.pi()]
            stats_pop["load"] += [mutation_load(fs_i, s, h)]
            stats_pop["efficacy"] += [efficacy_of_selection(fs_i, s, h)]

    return pandas.DataFrame.from_dict(stats_pop)


def mutation_load(fs, s, h):
    mu_1 = onepop_moments(fs, 1)
    if(isclose(h, 0.5)): 
        return s * mu_1 
    else:
        mu_2 = onepop_moments(fs, 2)
        return s * (2 * h * mu_1 + (1 - 2 * h) * mu_2)

def pi_k(fs, k):
    return 2. * (onepop_moments(fs, k) - onepop_moments(fs, k-1))

def Gamma_kh(fs, k, h):
    if(isclose(h, 0.5)): 
        return pi_k(fs, k) 
    else:
        return 2. * (h * pi_k(fs, k) + (1. - 2. * h) * pi_k(fs, k + 1))

def efficacy_of_selection(fs, s, h):
    if(isclose(h, 0.5)): 
        return s * s * Gamma_kh(fs, 1, h) / 4.
    else:
        w =  0.5 * h * Gamma_kh(fs, 1, h) + (1. - 2. * h) * Gamma_kh(fs, 2, h)
        return  s * s * w


mutation_load(fs, 0.1, 0.5)
efficacy_of_selection(fs, 0.1, 0.3)
pi_k(fs, 2)

params = np.array([2.10065897, 0.25066579, 0.22247642, 3.05297944,
                   0.09022469, 5.82773903, 3.79104318, 0.25730946,
                   0.12569788, 1.07182332, 0.36429414, 0.1108222, 
                   0.07072507])

(nuAf, nuB, nuEu0, nuEu, nuAs0, nuAs, 
 mAfB, mAfEu, mAfAs, mEuAs, TAf, TB, TEuAs) = params

fs_history = OutOfAfrica_stepwise(params, (4, 4, 4), 20)

fs_history[1][1].pi()

#fs_admix = moments.Manips.admix_into_new(fs, 0, 1, 2, 0.5)

#sts = moments.LinearSystem_1D.steady_state_1D(9, 1)
#fs = moments.Spectrum(sts)
#fs = moments.Manips.split_1D_to_2D(fs, 3, 6)
#fs = moments.Manips.split_2D_to_3D_2(fs, 3, 3)
