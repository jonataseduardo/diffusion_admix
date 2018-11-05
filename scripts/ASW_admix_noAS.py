import moments
import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
import OutOfAfrica as ooa
import matplotlib.pyplot as plt
import argparse
from gen_load_ultis import *


def ASW_admix_noAS((nuAf, nuB, nuEu0, nuEu,  
                    mAfB, mAfEu, TAf, TB, TEuAs), 
                   (n1, n2, n3),
                   (gamma1, gamma2),
                   (h1, h2), 
                   fadmix
                  ):
    """
    A three-population model used to model out-of-Africa demography.
    """
    #first step: a single population
    sts = moments.LinearSystem_1D.steady_state_1D(n1 + n2 + 2 * n3)
    fs = moments.Spectrum(sts)
    #integrate for time TAf (with constant population)
    fs.integrate([nuAf], TAf, 0.05, gamma = [gamma1], h = [h1])
    
    #separate into two populations.
    fs = moments.Manips.split_1D_to_2D(fs, n1 + n3, n2 + n3)
    
    #integrate two populations
    # migration rates matrix
    mig1=numpy.array([[0, mAfB],[mAfB, 0]])
    fs.integrate([nuAf, nuB], TB, 0.05, 
                 m=mig1, gamma = [gamma1, gamma2], h = [h1, h2])
    
    #define functions for population sizes
    # Removing Asian population to speed up simulations
    nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/TEuAs)
    nu2 = lambda t: [nuAf, nuEu_func(t)]

    # migration rates matrix
    fs.integrate(nu2, TEuAs, 0.05, m=mig1, 
                 gamma = [gamma1, gamma2], h = [h1, h2])
                                
    fs = moments.Manips.admix_into_new(fs, 0, 1, n3, fadmix)

    return fs


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "African Amiericans Admixture with moments")

    parser.add_argument('-ns', '--nsamples', 
                        type = int, nargs = '?', const = 1, default = 10, 
                        help = 'Number of samples for each population.') 

    parser.add_argument('-g', '--gamma',  
                        type = float, nargs = '?', const = 1, default = 0, 
                        help = 'Selection coef scaled by the initial population size (N_0 s).') 

    parser.add_argument('-d', '--dominance', 
                        type = float, nargs = '?', const = 1, default = 0.5, 
                        help = 'Dominance coeficient (0 = recessive, 0.5 = additive, 1 = dominance).') 

    parser.add_argument('-f', '--fadmix', 
                        type = float, nargs = '?', const = 1, default = 0.759, 
                        help = 'Proportion of the African ancestry') 

    parser.add_argument('-o', '--out', 
                        type = str, nargs = '?', const = 1, default = None, 
                        help = 'Output file path/filename_prefix') 

    args = parser.parse_args()

    nsamples = args.nsamples
    gamma = args.gamma
    h = args.dominance
    fadmix = args.fadmix

    (n1, n2, n3) = [nsamples] * 3
    (gamma1, gamma2) = [gamma] * 2
    (h1, h2) = [h] * 2
    
    params = np.array([2.10065897, 0.25066579, 0.22247642, 3.05297944,
                       3.79104318, 0.25730946, 0.36429414, 0.1108222, 
                       0.07072507])

    (nuAf, nuB, nuEu0, nuEu, 
     mAfB, mAfEu, TAf, TB,
     TEuAs) = params

    fs = ASW_admix_noAS((nuAf, nuB, nuEu0, nuEu, 
                         mAfB, mAfEu,TAf, TB, TEuAs), 
                        (n1, n2, n3),
                        (gamma1, gamma2),
                        (h1, h2), 
                        fadmix
                        )

    if(args.out is None):
        fn = 'ASW_ns{}_g{}_d{}_f{}.fs'.format(nsamples, gamma, h, fadmix)
        fs.to_file(fn)
    else:
        fs.to_file(args.out)

