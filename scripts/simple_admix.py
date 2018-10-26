import moments
import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
import OutOfAfrica as ooa
import seaborn as sea
import matplotlib.pyplot as plt
from gen_load_ultis import *


s = [0.5, 0.5, 0.5]
h = [0.5, 0.5, 0.5]
n = (4, 4, 4)
steps =(5, 5, 10) 

%time fs_history = OutOfAfrica_stepwise(s, h, n, steps)

fs_stats = OOA_stats(fs_history, s, h) 
fs_stats
fs_stats[fs_stats.time > 0.5]

sea.set()

sea.relplot(data = fs_stats[fs_stats.time > 0.5],
            x = 'time', 
            y = 'load', 
            kind = 'line',
            hue = 'pop_id')

plt.show()

#fs_admix = moments.Manips.admix_into_new(fs, 0, 1, 2, 0.5)

sts = moments.LinearSystem_1D.steady_state_1D(9, 1)
fs = moments.Spectrum(sts)
# cross_moments(fs,1)
# fs = moments.Manips.split_1D_to_2D(fs, 3, 6)
# fs = moments.Manips.split_2D_to_3D_2(fs, 3, 3)
