import slim_output_parser as sop
import pandas as pd
import matplotlib.pyplot as plt
import os
import moments 
import seaborn as sns

def mutation_load(m_dt):
    h = m_dt["dominance"]
    s = m_dt["selection"]
    return s * (2. * h * m_dt["mu_1"] + (1. - 2. * h) * m_dt["mu_2"])

def pi_k(m_dt, k):
    return 2. * (m_dt["mu_" + str(k)] - m_dt["mu_" + str(k + 1)])

def htz(m_dt):
    return 2. * (m_dt["mu_2"] - m_dt["mu_1"] ** 2)

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
        m_dt["mu_" + str(k)] = eval_moment(m_dt, k)

    m_dt["fit"] = fit_efficacy(m_dt) 
    m_dt["morton"] = morton_efficacy(m_dt)
    m_dt["load"] = mutation_load(m_dt)
    m_dt["htz"] = htz(m_dt)

    return m_dt

def huber_dominance(m_dt, h_intercept = 0.5, h_rate = 1e6):
    s = m_dt["selection"] 
    s0 = 1. / h_intercept
    return 1. / (s0 - h_rate * s)

def get_simul_summaries(f, sel_bins = [-1, -0.01, -0.001, -0.0001, 0]):

    m_dt = sop.slim_output_parser(f)[0]

    if  "huber" in f:
        m_dt.dominance = huber_dominance(m_dt)
        m_dt['dominance_key'] = "huber"
    else:
        m_dt['dominance_key'] = m_dt.dominance.apply(str)

    eval_stats_from_slim(m_dt, inplace = True)

    intervalindex = sel_bins
    m_dt["selection_bin"] = pd.cut(m_dt.selection, 
                                    intervalindex)
                                    #labels = False)

    stats = {"selection": ["count", "mean"], 
             "htz": ["sum", "mean"],
             "load": ["sum", "mean"],
             "fit": ["sum", "mean"],
             "morton": ["sum", "mean"],
             "mu_1": ["sum", "mean"],
             "mu_2": ["sum", "mean"],
             "mu_3": ["mean", "sum"],
             }

    run_s = m_dt.groupby(["dominance_key", "focal_pop_id", "selection_bin"]
                        ).agg(stats)

    run_s.columns = ['_'.join(col) for col in run_s.columns]

    rid = get_rid(f)
    run_s["rid"] = rid
    return run_s.reset_index()

def eval_proportion_to_pop(df_stats, pop = "p1"):
    df = df_stats.sort_values(by = ["level_0", "focal_pop_id", "selection_bin"])
    df.index = range(df.shape[0])

    d0 = df[df.focal_pop_id == pop]
    d0.index = range(d0.shape[0])
    dd0 = pd.concat([d0] * df.focal_pop_id.unique().shape[0])
    dd0 = dd0.sort_values(by = ["level_0", "focal_pop_id", "selection_bin"])
    dd0.index = range(dd0.shape[0])


    ratios = df.iloc[:,3:] / dd0.iloc[:,3:]
    ratios = pd.concat([df.iloc[:,:3],
                        ratios], 
                       axis = 1)

    return ratios

def get_rid(f):
    return int(f.split("-")[-1].split('.')[0])

def dominance_ticks(ss):
    ss.loc[ss.dominance_key == "huber", "dominance_key"] = "-0.1"
    ss.dominance_key = pd.to_numeric(ss.dominance_key)
    xticks = [str(i) for i in sorted(ss.dominance_key.unique())]
    xticks[0] = 'huber'
    return xticks

def boxplot_slim_bysel(slim_summaries,
                       summ_col = 'mu_1_sum',
                       ylabel = None,
                       show = True):

    ss = slim_summaries.copy()
    ticks = dominance_ticks(ss)

    ax = sns.catplot(data = ss,
                     kind = 'box',
                     col = 'selection_bin',
                     col_wrap = 2,
                     sharey = False,
                     x = 'dominance_key',
                     y = summ_col,
                     hue = "focal_pop_id")

    ax.set_xlabels("dominance")
    ax.set_xticklabels(ticks)
    if ylabel is not None:
        ax.set_ylabels(ylabel)

    fname = summ_col + '_bs.png'
    ax.savefig('../figures/' + fname)
    if show:
        plt.show()

def lineplot_slim_bysel(slim_summaries,
                        summ_col = 'mu_1_sum',
                        stat = 'mean',
                        ylabel = None,
                        show = False):

    stat_s = {summ_col : stat}
    ss = slim_summaries.groupby(['selection_bin', 'dominance_key', 'focal_pop_id']
                                ).agg(stat_s).reset_index()
    ticks = dominance_ticks(ss)

    #sns.set(style="white")

    #palette = dict(zip(ratios.fadmix.unique(),
    #                   sns.color_palette("rocket_r", len(ratios.fadmix.unique()))))

    ax = sns.relplot(data = ss, 
                     x = "dominance_key", 
                     y = summ_col,
                     hue = "focal_pop_id",
                     col = 'selection_bin',
                     col_wrap = 2,
                     kind = "line",
                     facet_kws = {"sharey" : False})

    ax.set_xlabels("dominance")
    #ax.set_xticklabels(ticks)

    ax.savefig('../figures/' + summ_col + '_line.png')

    if show: 
        plt.show()


def boxplot_slim_full(slim_summaries, 
                      summ_col = 'mu_1_sum', 
                      stat = 'sum',
                      ylabel = None,
                      show = False):

    stat_s = {summ_col : stat}
    ss = slim_summaries.groupby(['rid', 'dominance_key', 'focal_pop_id']
                                ).agg(stat_s).reset_index()
    ticks = dominance_ticks(ss)

    ax = sns.boxplot(data = ss, 
                     x = 'dominance_key', 
                     y = summ_col,
                     hue = "focal_pop_id") 

    ax.set(xlabel = "dominance")
    ax.set_xticklabels(ticks)

    if ylabel is not None:
        ax.set(ylabel = ylabel)

    fig = ax.get_figure()
    if show:
        fig.show()

    fname = summ_col + '_full.png'
    fig.savefig('../figures/' + fname)



if __name__ == "__main__":


    path = '../slim_data/'
    list_files = os.listdir(path)

    fns = [f for f in list_files if 'l-1000000_' in f]
    fh = [f for f in fns if ('huber' not in f) and ('neutral' not in f) ]

    summaries = pd.concat([get_simul_summaries(f) for f in fh])

    stats_s = {'load_sum': 'sum', 'mu_1_sum': 'sum', 'selection_count': 'sum'}
    s = summaries.groupby(['rid', 'dominance_key','focal_pop_id']
            ).agg(stats_s).reset_index()


    summaries.groupby(["dominance", "focal_pop_id"]).load_sum.sum()
    summaries.groupby(["dominance", "focal_pop_id"]).selection_count.mean()
    summaries.groupby(["dominance", "focal_pop_id"]).mu_1_sum.sum()

    summaries.groupby("focal_pop_id").mu_1_sum.sum()
    
    summaries[summaries.dominance_key == 'huber'] 
    ss = s[(s.dominance_key == "huber")]
    ss = s[(s.focal_pop_id == "p1")]

    ax = sns.boxplot(data = ss, 
                       x = 'dominance_key', 
                       y = 'mu_1_sum') 

    plt.show()
