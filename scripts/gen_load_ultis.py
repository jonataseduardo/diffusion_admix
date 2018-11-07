import numpy

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def onepop_moments(fs, k):
    ns = fs.sample_sizes[0]
    # sample frequencies p 
    p = numpy.arange(0, ns + 1, dtype=float) / ns
    return numpy.ma.sum(fs *(p ** k))

def cross_moments(fs, k):
    r = fs.Npop
    ns = fs.sample_sizes

    # counts_per_pop is an r+1 dimensional array, where the last axis simply
    # records the indices of the entry. 
    # For example, counts_per_pop[4,19,8] = [4,19,8]
    counts_per_pop = numpy.indices(fs.shape)
    counts_per_pop = numpy.transpose(counts_per_pop, axes=range(1, r + 1) + [0])
    
    ktilde = numpy.tile(k, numpy.append(ns + 1, 1))

    x = 1. * counts_per_pop / ns 

    x_k = numpy.power(x, ktilde)
    
    this_slice = tuple([slice(None)] * r + [numpy.newaxis])

    fstilde = numpy.tile(fs.data[this_slice], [1] * r  + [r])

    mu_k = x_k * fstilde

    return numpy.na.sum(mu_k)

def mutation_load(fs, s, h):
    mu_1 = onepop_moments(fs, 1)
    if(isclose(h, 0.5)): 
        return s * mu_1 
    else:
        mu_2 = onepop_moments(fs, 2)
        return s * (2 * h * mu_1 + (1 - 2 * h) * mu_2)

def pi_k(fs, k):
    return 2. * (onepop_moments(fs, k) - onepop_moments(fs, k + 1))

def Gamma_kh(fs, k, h):
    if(isclose(h, 0.5)): 
        return pi_k(fs, k) 
    else:
        return 2. * (h * pi_k(fs, k) + (1. - 2. * h) * pi_k(fs, k + 1))

def fit_efficacy(fs, s, h):
    if(isclose(h, 0.5)): 
        return s * s * Gamma_kh(fs, 1, h) / 4.
    else:
        w =  0.5 * h * Gamma_kh(fs, 1, h) + (1. - 2. * h) * Gamma_kh(fs, 2, h)
        return  s * s * w

def morton_efficacy(fs, s, h):
    return i * s * Gamma_kh(fs, 1, h) / 4.
