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

def huber_dominance(m_dt, h_intercept = 0.5, h_rate = 1e6):
    s = m_dt["selection"] 
    s0 = 1. / h_intercept
    return s / (s0 + h_rate * s)

def get_simul_summaries(f):

    m_dt = sop.slim_output_parser(path + f)[0]

    if  "huber" in f:
        m_dt.dominance = huber_dominance(m_dt)

    eval_stats_from_slim(m_dt, inplace = True)

    intervalindex = [-1.0, -0.1, -0.01, -0.001, -0.0001, 0.0]
    m_dt["selection_bins"] = pd.cut(m_dt.selection, 
                                    intervalindex,
                                    labels = False)

    run_summary = m_dt.groupby(["focal_pop_id", "selection_bins"]
                                ).agg({"selection": ["count", "mean"], 
                                    "load_i": ["sum", "mean"],
                                    "fit_i": ["sum", "mean"],
                                    "morton_i": ["sum", "mean"],
                                    "mu_i1": ["sum", "mean"],
                                    "mu_i2": ["sum", "mean"],
                                    "mu_i3": ["mean", "sum"],
                                    })

    return run_summary

if __name__ == "__main__":


    path = '../data/'
    list_files = os.listdir(path)

    N0 = 11000

    fns = [f for f in list_files if 'huber_l-1000000_' in f]

    %time summaries = [ get_simul_summaries(f) for f in fns]

    summaries[1].load_i
    s = pd.concat(summaries, keys = range(len(fns)))
    s
    s.columns
    s.load_i


