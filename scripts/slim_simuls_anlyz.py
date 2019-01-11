import slim_output_parser as sop
import pandas as pd
import matplotlib.pyplot as plt
import os
import moments 
import admix_and_analyse as aaa


def mutation_load(m_dt):
    h = m_dt["dominance"]
    s = m_dt["selection"]
    return s * (2. * h * m_dt["mu_i1"] + (1. - 2. * h) * m_dt["mu_i2"])

def pi_k(m_dt, k):
    return 2. * (m_dt["mu_i" + str(k)] - m_dt["mu_i" + str(k + 1)])

def Gamma_kh(m_dt, k):
    h = m_dt["dominance"]
    return 2. * (h * pi_k(m_dt, k) + (1. - 2. * h) * pi_k(m_dt, k + 1))

def fit_efficacy(m_dt):
    h = m_dt["dominance"]
    s = m_dt["selection"]
    w =  0.5 * h * Gamma_kh(m_dt, 1) + (1. - 2. * h) * Gamma_kh(m_dt, 2)
    return  s * s * w

def morton_efficacy(m_dt):
    s = m_dt["selection"]
    return 0.25 * s * s * Gamma_kh(m_dt, 1) 

def eval_moment(m_dt, k):
    return m_dt.allele_frequency ** k

def eval_stats_from_slim(mutation_db, inplace = True):
    if inplace:
        m_dt = mutation_db
    else:
        m_dt = mutation_db.copy()

    m_dt["allele_frequency"] = m_dt.allele_count / m_dt.sample_size
    for k in range(1,6): 
        m_dt["mu_i" + str(k)] = eval_moment(m_dt, k)

    m_dt["fit_i"] = fit_efficacy(m_dt) 
    m_dt["morton_i"] = morton_efficacy(m_dt)
    m_dt["load_i"] = mutation_load(m_dt)

    return m_dt


path = '../data/'
list_files = os.listdir(path)

N0 = 11000

fns = [f for f in list_files if 'huber_l-1000000_' in f]

f = fns[35]

m_dt = sop.slim_output_parser(path + f)[0]
m_dt = eval_stats_from_slim(m_dt, inplace = True)

m_dt["selection_bins"] = pd.cut(m_dt.selection, 4, labels = False)

m_dt.iloc[1:40,[2,3,4,5,6,8,9, -2, -1]]

m_dt.groupby(["focal_pop_id", "selection_bins"]).selection.mean()
m_dt.groupby(["focal_pop_id", "selection_bins"]
            ).agg({"selection": "mean", 
                   "load_i": "mean",
                   "fit_i": "mean",
                   "morton_i": "mean",
                   "mu_i1": "mean",
                   "mu_i2": "mean",
                   "mu_i3": "mean"
                   })

ss = sop.make_slim_sfs(m_dt_cut)[0]
ms = moments.Spectrum(ss)
mss = ms.project([200])

aaa.stats_from_fs(mss, 1., 0.5)

s = mss.marginalize([1])

s.pi()
