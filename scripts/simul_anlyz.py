import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import moments as mm 
from itertools import combinations
import gen_load_ultis as gl
import os, sys

def get_fname_parameters(fname):
    fp = fname[:-3].split("_")[1:]
    g = h = f = None
    for p in fp:
        if p.startswith('g'):
            g = float(p[1:])
        if p.startswith('d'):
            h = float(p[1:])
        if p.startswith('f'):
            f = float(p[1:])
    return {"gamma": g, "dominance": h, "fadmix": f}

def stats_from_fs(fs, gamma, h):
    """
    Read Simulation SFS file and eval the k first moments for of each
    population.

    Params: 
        fs  Spectrum object
        gamma  selection parameter
        h Dominance
    Return:
        Pandas.DataFrame with multple statistics 
    """

    stats = dict()
    stats["pop_id"] = []
    stats["pi"] = []
    stats["load"] = []
    stats["fit"] = []
    stats["morton"] = []
    stats["gamma"] = gamma
    stats["dominance"] = h

    for k in range(5):
        stats["mu_" + str(k +1)] = []

    np = fs.Npop
    rnp = range(np)
    lpop = [(list(set(rnp) ^ set(i)), i) for i in combinations(rnp, np - 1)]
    for i in lpop:
        pid = i[0][0]
        stats["pop_id"] += [pid]
        try:
            fs_i = fs.marginalize(i[1])
        except:
            fs_i = fs

        stats["pi"] += [fs_i.pi()]
        stats["load"] += [gl.mutation_load(fs_i, gamma[pid], h[pid])]
        stats["fit"] += [gl.efficacy_of_selection(fs_i, gamma[pid], h[pid])]
        stats["morton"] += [gl.efficacy_of_selection(fs_i, gamma[pid], 0.5)]
        for k in range(5):
            stats["mu_" + str(k +1)] += [gl.onepop_moments(fs_i, k + 1)]
    
    dt_stats = pd.DataFrame(stats)
    
    c_order = ["pop_id", "gamma", "dominance", "pi", "load", "fit", "morton"]
    c_order += ["mu_" + str(k +1) for k in range(5)]
    return dt_stats[c_order]

def get_simul_stats(fpath, nsimuls = None):
    """
    Read Simulation SFS file and eval the k first moments for of each
    population.

    Params: 
        fpath file path: string
        k  the total number of moments to be evalueted
    Return:
        Pandas.DataFrame with multple statistics 
    """
    df_list = []

    list_files =  [f for f in os.listdir(fpath) if f.endswith(".fs")]

    if nsimuls is not None:
        list_files = list_files[:nsimuls] 

    for fn in list_files:
        fs = mm.Spectrum.from_file(fpath + fn)
        f_par = get_fname_parameters(fn)
        np = fs.Npop
        g = [f_par["gamma"]] * np
        h = [f_par["dominance"]] * np
        df = stats_from_fs(fs, g, h) 
        df["fadmix"] = [f_par["fadmix"]] * np
        df_list += [df]

    df_stats = pd.concat(df_list, ignore_index = True)

    return df_stats

if __name__ == "__main__":


    df_stats = get_simul_stats('../data/')

    df = df_stats.dropna().sort_values(by = ["pop_id", "gamma", 
                                             "dominance", "fadmix"])

    pop_dt = dict()
    for i in range(3):
        pop_dt[i]  = df[df.pop_id == i]
        pop_dt[i].index = range(pop_dt[i].shape[0]) 

    pop_dt[1].iloc[0:5,1:]

