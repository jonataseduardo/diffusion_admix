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
    stats["num_sites"] = [fs.S()]
    stats["load"] = [gl.mutation_load(fs, gamma, h)]
    stats["fit"] = [gl.fit_efficacy(fs, gamma, h)]
    stats["morton"] = [gl.morton_efficacy(fs, gamma, h)]
    stats["gamma"] = gamma
    stats["dominance"] = h
    for k in range(5):
        stats["mu_" + str(k +1)] = [gl.onepop_moments(fs, k + 1)]

    dt_stats = pd.DataFrame(stats)

    c_order = ["gamma", "dominance", "num_sites", "pi", "load", "fit", "morton"]
    c_order += ["mu_" + str(k +1) for k in range(5)]
    return dt_stats[c_order]

def admix_and_get_simul_stats(fpath, nsimuls = None, n = 3, file_pattern = None):
    """
    Read Simulation SFS file and eval the k first moments for of each
    population.

    Params: 
        fpath file path: string
        k  the total number of moments to be evalueted
        n number of partitions 
    Return:
        Pandas.DataFrame with multple statistics 
    """

    list_files =  sorted([f for f in os.listdir(fpath) if f.endswith(".fs")])

    if nsimuls is not None:
        list_files = list_files[:nsimuls] 

    if file_pattern is not None:
        list_files = [f for f in list_files if file_pattern in f]

    l_ratios = numpy.round(numpy.linspace(0, 1, n + 2), decimals = 2)[1:-1]
    df_list = []
    for fn in list_files:
        fs = moments.Spectrum.from_file(fpath + fn)
        # For extreme values of gamma numerical erros will return nans and
        # negative values in the sfs
        #fs[numpy.isnan(fs)] = 0 
        #fs[fs < 0] = 0 
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
        for f, pid in zip(l_ratios, range(2, len(l_ratios) + 2)):
            fs_admix = moments.Manips.admix_into_new(fs, 0, 1, admix_size, f)
            df_admix = stats_from_fs(fs_admix, g, h) 
            df_admix["pop_id"] = pid 
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

def eur_matrix_plot(ratios, col_name = "pi"):
    ratios = ratios[ratios.gamma < - 0.5]
    yticks = ratios.gamma.unique()
    xticks = ratios.dominance.unique()
    ncols = xticks.shape[0]
    nrows = yticks.shape[0]
    stats = ratios[ratios.pop_id == 1].loc[:, col_name]
    mstats = stats.values.reshape(nrows, ncols)

    smin =  stats.min()
    smax = stats.max() 
    dc = max(smax - 1, 1 - smin)
    max_color = 1 + dc
    min_color = 1 - dc

    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = ax.imshow(mstats, 
                   interpolation = 'spline16', 
                   cmap = 'RdBu',
                   vmin = min_color, 
                   vmax = max_color)
    ax.set_xticks(range(ncols))
    ax.set_xticklabels(xticks)
    ax.set_yticks(range(nrows))
    ax.set_yticklabels(yticks)
    cbar = fig.colorbar(im)
    cbar.ax.set_ylabel(col_name)
    fig.tight_layout()
    fig.savefig("../figures/pdiag_" + col_name + ".png")
    return fig

def plot_admix_hueline(ratios, stats_key): 

    sns.set(style="white")


    palette = dict(zip(ratios.fadmix.unique(),
                       sns.color_palette("rocket_r", len(ratios.fadmix.unique()))))

    sns_plot = sns.relplot(data = ratios[ratios.gamma < -0.2], 
                           x = "dominance", 
                           y = stats_key,
                           hue = "fadmix",
                           col = "gamma",
                           col_wrap = 3,
                           palette = palette, 
                           kind = "line")
    sns_plot.savefig('../figures/' + stats_key + '.png')

if __name__ == "__main__":

    df_stats = admix_and_get_simul_stats('../moments_data/')

    ratios = eval_proportion_to_Afr(df_stats)
    
    plot_admix_hueline(ratios, "mu_1")
    plot_admix_hueline(ratios, "mu_2")
    plot_admix_hueline(ratios, "pi")
    plot_admix_hueline(ratios, "fit")
    plot_admix_hueline(ratios, "morton")
    plot_admix_hueline(ratios, "load")

    eur_matrix_plot(ratios, "mu_1")
    eur_matrix_plot(ratios, "mu_2")
    eur_matrix_plot(ratios, "pi")
    eur_matrix_plot(ratios, "fit")
    eur_matrix_plot(ratios, "morton")
    eur_matrix_plot(ratios, "load")

    fig = eur_matrix_plot(ratios[ratios.pop_id == 1], "load")
    fig.show()
