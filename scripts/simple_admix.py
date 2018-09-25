import moments
import matplotlib.pyplot as pl
import numpy as np
import numpy
import OutOfAfrica as ooa


params = np.array([2.10065897, 0.25066579, 0.22247642, 3.05297944,
                   0.09022469, 5.82773903, 3.79104318, 0.25730946,
                   0.12569788, 1.07182332, 0.36429414, 0.1108222, 
                   0.07072507])


def test(n, theta):
    sts = moments.LinearSystem_1D.steady_state_1D(n, theta = theta)
    fs = moments.Spectrum(sts)

    def mu_k_1d(fs, k):
        ns = fs.sample_sizes[0]
        # sample frequencies p 
        p = np.arange(0, ns + 1, dtype=float) / ns
        return (fs *(p ** k)).sum()

    print(n, mu_k_1d(fs, 1), mu_k_1d(fs, 2))


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
    
    this_slice = [slice(None)] * r + [numpy.newaxis]

    fstilde = numpy.tile(fs[this_slice], [1] * r  + [r])

    mu_k = x_k * fstilde

    return mu_k.sum()



# sts = moments.LinearSystem_1D.steady_state_1D(9, 1)
# fs = moments.Spectrum(sts)
# fs = moments.Manips.split_1D_to_2D(fs, 3, 6)
# fs = moments.Manips.split_2D_to_3D_2(fs, 3, 3)



(n1, n2, n3) = (4, 4, 4)

(nuAf, nuB, nuEu0, nuEu, nuAs0, nuAs,mAfB, mAfEu, mAfAs, mEuAs, TAf, TB, TEuAs) = params
fs = ooa.OutOfAfrica(params, (4, 4, 4))

fs_admix = moments.Manips.admix_into_new(fs, 0, 1, 2, 0.5)
