import slim_output_parser as sop
import pandas as pd
import matplotlib.pyplot as plt
import os
import moments 

path = '../data/'
list_files = os.listdir(path)

fns = [f for f in list_files if 'huber_l-1000000_' in f]

f = fns[35]

m_dt = sop.slim_output_parser(path + f)[0]

m_dt["selection_bins"] = pd.cut(m_dt.selection, 4, labels = False)

m_dt_cut = m_dt[(m_dt["selection_bins"] == 2) & 
                (m_dt["focal_pop_id"].isin(["p1", "p2"]))]

m_dt.groupby(["selection_bins"]).selection.median()

ss = sop.make_slim_sfs(m_dt_cut)[0]
ms = moments.Spectrum(ss)
mss = ms.project([200,200])

s = mss.marginalize([1])

s.pi()
