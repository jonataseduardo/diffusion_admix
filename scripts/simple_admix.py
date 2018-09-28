import moments
import matplotlib.pyplot as pl
import numpy as np
import numpy
import OutOfAfrica as ooa
from collections import OrderedDict

def onepop_moments(fs, k):
    ns = fs.sample_sizes[0]
    # sample frequencies p 
    p = np.arange(0, ns + 1, dtype=float) / ns
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

    time_sfs = OrderedDict()

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
            time_sfs[step] = [time, fs]
            step += 1

        if(time > TAf and time <= 1.001 * (TAf + TB)):
            if(split_1):
                print('first split' , step, time )
                split_1 = False
                dt = dtB
                fs = moments.Manips.split_1D_to_2D(fs, n1, n2+n3)
                time_sfs[step] = [time, fs]
                step += 1
                mig1=numpy.array([[0, mAfB],[mAfB, 0]])

            print('after first split' , time )
            fs.integrate([nuAf, nuB], dt, 0.5 * dt, m=mig1)
            time_sfs[step] = [time, fs]
            step += 1

        if(time > TAf + TB and time <= 1.001 * (TAf + TB + TEuAs)):
            if(split_2):
                print('sencond split' , step, time )
                split_2 = False
                dt = dtEuAs
                fs = moments.Manips.split_2D_to_3D_2(fs, n2, n3)
                time_sfs[step] = [time, fs]
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

            time_sfs[step] = [time, fs]
            nuEu0t = nuEu0t*(nuEu/nuEu0t)**(dt/TEuAs)
            nuAs0t = nuAs0t*(nuAs/nuAs0t)**(dt/TEuAs)
            step += 1
                                
    return time_sfs

params = np.array([2.10065897, 0.25066579, 0.22247642, 3.05297944,
                   0.09022469, 5.82773903, 3.79104318, 0.25730946,
                   0.12569788, 1.07182332, 0.36429414, 0.1108222, 
                   0.07072507])

(nuAf, nuB, nuEu0, nuEu, nuAs0, nuAs,mAfB, mAfEu, mAfAs, mEuAs, TAf, TB, TEuAs) = params
fs_t = OutOfAfrica_stepwise(params, (4, 4, 4), 20)

(5. / (5. - 1.)) * 2 * (cross_moments(fs_t[1][1], 1) - cross_moments(fs_t[1][1], 2))
fs_t[1][1].pi()

fs_admix = moments.Manips.admix_into_new(fs, 0, 1, 2, 0.5)

sts = moments.LinearSystem_1D.steady_state_1D(9, 1)
fs = moments.Spectrum(sts)
n = fs.sample_sizes
(n / (n - 1.)) * 2. * (cross_moments(fs, 1) - cross_moments(fs, 2))
fs.pi()

#fs = moments.Manips.split_1D_to_2D(fs, 3, 6)
#fs = moments.Manips.split_2D_to_3D_2(fs, 3, 3)
