import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import slim_simul_anlz as ssa
import moments_simul_anlz as msa
import os
import moments
from importlib import reload
from scipy import interpolate
import scipy

reload(msa)
fpath = "../moments_data/"

ff = [f for f in list_files if "g-1000" in f]
fn = ff[0]
fn

        fs = moments.Spectrum.from_file(fpath + fn)
fs[np.isnan(fs)] = 0 
fs[fs < 0] = 0 
fs

def gamma_dist(mgamma, alpha, beta):
    """
    x, shape, scale
    """
    return scipy.stats.distributions.gamma.pdf(-mgamma, alpha, scale=beta)    

def gamma_cdf(mgamma, alpha, beta):
    """
    x, shape, scale
    """
    return scipy.stats.distributions.gamma.cdf(-mgamma, alpha, scale=beta)    


if __name__ == "__main__":

    #Slim3 simulation parameters
    N_0 = 11273
    alpha = -0.00657402 
    beta = 0.186
    
    slim_files = os.listdir('../slim_data/')

    #slim files list
    sfl = [ f for f in slim_files if ('l-1000000_' in f) and ('neutral' not in f)]

    slim_summaries = pd.concat([ssa.get_simul_summaries('../slim_data/' +f) for f in sfl])

    stats_s = {'load_sum': 'sum', 'mu_1_sum': 'sum', 'selection_count': 'sum'}
    s = slim_summaries.groupby(['rid', 'dominance_key','focal_pop_id']
                               ).agg(stats_s).reset_index()

    slim_summaries.loc[slim_summaries.dominance_key == "huber", "dominance_key"] = "-0.05"
    slim_summaries["dominance_key"] = slim_summaries.dominance_key.apply(np.float)

    slim_summaries.groupby(["dominance_key", "focal_pop_id"]).load_sum.sum()
    slim_summaries.groupby(["dominance_key", "focal_pop_id"]).selection_count.mean()
    slim_summaries.groupby(["dominance_key", "focal_pop_id"]).mu_1_sum.sum()

    N_0 * 2 * slim_summaries[slim_summaries.selection_bins == 4].selection_mean.mean()
    
    ss = s[(s.dominance_key == "huber")]
    ss = s[(s.focal_pop_id == "p1")]

    ss.dominance_key.unique()

    ax = sns.boxplot(data = ss, 
                       x = 'dominance_key', 
                       y = 'mu_1_sum') 

    plt.show()

    len(df_stats.gamma.unique())

    z = df_stats.loc[(df_stats.pop_id == 0) & (df_stats.dominance == 0.0)]
    x = z.gamma
    y = z.pi

    f = interpolate.interp1d(x, y, kind = 'linear')
    xnew = np.linspace(min(x), max(x), 50)
    ynew = f(xnew)
    plt.plot(x,y, 'o', xnew, ynew, '-')
    plt.show()

    reload(msa)
    df_stats = msa.admix_and_get_simul_stats('../moments_data/')

    ratios = msa.eval_proportion_to_Afr(df_stats)
    ratios[np.isnan(ratios)] = 1

    ratios[ratios.gamma == -1000]
    df_stats[df_stats.gamma == -1000].iloc[:,1:8]
    df_stats1[df_stats1.gamma == -1000].iloc[:,1:8]
    df_stats1 = df_stats.copy()
    df_stats[df_stats.abs() < 1e-10] = 0
    
    df_stats

    fig = msa.eur_matrix_plot(ratios[ratios.pop_id == 1], "load")
    msa.plot_admix_hueline(ratios[(ratios.gamma > -10)],
                                  #(ratios.dominance <0.3) & 
                                  #(ratios.gamma > -100)], 
                           "mu_1")
    
    fig.show()
