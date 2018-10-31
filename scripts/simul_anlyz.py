import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import moments as mm 
import moments 
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


def stats_from_file(fname, gamma, h):
    """
    Read Simulation SFS file and eval the k first moments for of each
    population.

    Params: 
        fname file name: string
        k  the total number of moments to be evalueted
    """
    fs = mm.Spectrum.from_file(fname)

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

if __name__ == "__main__":


    list_files = os.listdir("../data")
    list_files
