import pandas as pd
import numpy
import moments
import matplotlib.pyplot as plt
from itertools import combinations
import gen_load_ultis as gl
import os, sys
import seaborn as sns

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
    stats["pi"] = [fs.pi()]
    stats["load"] = [gl.mutation_load(fs, gamma, h)]
    stats["fit"] = [gl.efficacy_of_selection(fs, gamma, h)]
    stats["morton"] = [gl.efficacy_of_selection(fs, gamma, 0.5)]
    stats["gamma"] = gamma
    stats["dominance"] = h
    for k in range(5):
        stats["mu_" + str(k +1)] = [gl.onepop_moments(fs, k + 1)]

    dt_stats = pd.DataFrame(stats)

    c_order = ["gamma", "dominance", "pi", "load", "fit", "morton"]
    c_order += ["mu_" + str(k +1) for k in range(5)]
    return dt_stats[c_order]

def admix_and_get_simul_stats(fpath, nsimuls = None):
    """
    Read Simulation SFS file and eval the k first moments for of each
    population.

    Params: 
        fpath file path: string
        k  the total number of moments to be evalueted
    Return:
        Pandas.DataFrame with multple statistics 
    """

    list_files =  sorted([f for f in os.listdir(fpath) if f.endswith(".fs")])

    if nsimuls is not None:
        list_files = list_files[:nsimuls] 

    l_ratios = numpy.round(numpy.linspace(0.1, 0.9, 9), decimals = 1)
    df_list = []
    for fn in list_files:
        fs = moments.Spectrum.from_file(fpath + fn)
        f_par = get_fname_parameters(fn)
        g = f_par["gamma"]
        h = f_par["dominance"]
        df0 = stats_from_fs(fs.marginalize([1]), g, h) 
        df1 = stats_from_fs(fs.marginalize([0]), g, h) 
        df0["pop_id"] = 0 
        df1["pop_id"] = 1 
        df0["fadmix"] = numpy.round(1.0, decimals = 1)
        df1["fadmix"] = numpy.round(0.0, decimals = 1)
        df_list += [df1, df0]
        admix_size = int(0.5 * numpy.array(fs.shape).sum()) - 1
        for f in l_ratios:
            fs_admix = moments.Manips.admix_into_new(fs, 0, 1, admix_size, f)
            df_admix = stats_from_fs(fs_admix, g, h) 
            df_admix["pop_id"] = 2 
            df_admix["fadmix"] = f 
            df_list += [df_admix]

    df_stats = pd.concat(df_list, ignore_index = True)

    return df_stats

def eval_proportion_to_Afr(df_stats):
    df = df_stats.sort_values(by = ["gamma", "dominance", "fadmix"])
    df.index = range(df.shape[0])

    d0 = df[df.pop_id == 0]
    d0.index = range(d0.shape[0])
    dd0 = pd.concat([d0] * df.fadmix.unique().shape[0])
    dd0 = dd0.sort_values(by = ["gamma", "dominance", "fadmix"])
    dd0.index = range(dd0.shape[0])

    #ratios = df.iloc[:,:-2] / dd0.iloc[:,:-2]
    ratios = df.iloc[:,2:-2] / dd0.iloc[:,2:-2]
    ratios = pd.concat([df.loc[:,["gamma", "dominance", "fadmix", "pop_id"]],
                        ratios], 
                       axis = 1)
    return ratios

if __name__ == "__main__":


    #fpath = "../data/"
    df_stats = admix_and_get_simul_stats('../data/')

    ratios = eval_proportion_to_Afr(df_stats)


    sns.set(style="ticks")

    palette = dict(zip(ratios.fadmix.unique(),
		       sns.color_palette("rocket_r", len(ratios.fadmix.unique()))))

    sns.relplot(data = ratios[ratios.gamma > -2], 
                x = "dominance", 
                y = "morton",
                hue = "fadmix",
		col = "gamma",
		palette = palette, 
                kind = "line")
    plt.show()




