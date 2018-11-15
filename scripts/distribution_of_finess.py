import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sps
import scipy.stats as sts

# Parameters from Fitdadi paper
# Kim et al, Genetics 345-361 (2017)
Nanc = 7778.0
Ncur = 3.352 * Nanc
beta = 0.179
alpha = 3161 / ( 2. * Nanc )

powers = - np.arange(4, 0, -1)
s =  np.power([10.0] * len(powers), powers )
cdf = sts.gamma.cdf(s, a = alpha, scale = beta)

Kim_cdf = np.array([27.7, 17.1, 20.9, 37.3])

x = np.linspace(sts.gamma.ppf(0.01, a = alpha, scale = beta),
                sts.gamma.ppf(0.99, a = alpha, scale = beta),
                1000)

fig, ax = plt.subplots(1, 1)
ax.plot(x, sts.gamma.pdf(x, a = alpha, scale = beta))
ax.set_ylim((0, 10))
fig.show()
