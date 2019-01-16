import slim_output_parser as sop
import pandas as pd
import matplotlib.pyplot as plt
import os
import moments 
import admix_and_analyse as aaa
import seaborn as sns


def mutation_load(m_dt):
    h = m_dt["dominance"]
    s = m_dt["selection"]
    return s * (2. * h * m_dt["mu_1"] + (1. - 2. * h) * m_dt["mu_2"])

def pi_k(m_dt, k):
    return 2. * (m_dt["mu_" + str(k)] - m_dt["mu_" + str(k + 1)])

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

    return m_dt

def huber_dominance(m_dt, h_intercept = 0.5, h_rate = 1e6):
    s = m_dt["selection"] 
    s0 = 1. / h_intercept
    return 1. / (s0 - h_rate * s)

def get_simul_summaries(f):

    m_dt = sop.slim_output_parser(path + f)[0]

    if  "huber" in f:
        m_dt.dominance = huber_dominance(m_dt)
        m_dt['dominance_key'] = "huber"
    else:
        m_dt['dominance_key'] = m_dt.dominance.apply(str)

    eval_stats_from_slim(m_dt, inplace = True)

    intervalindex = [-1.0, -0.1, -0.01, -0.001, -0.0001, 0.0]
    m_dt["selection_bins"] = pd.cut(m_dt.selection, 
                                    intervalindex,
                                    labels = False)

    stats = {"selection": ["count", "mean"], 
             "load": ["sum", "mean"],
             "fit": ["sum", "mean"],
             "morton": ["sum", "mean"],
             "mu_1": ["sum", "mean"],
             "mu_2": ["sum", "mean"],
             "mu_3": ["mean", "sum"],
             }

    run_s = m_dt.groupby(["dominance_key", "focal_pop_id", "selection_bins"]
                        ).agg(stats)

    run_s.columns = ['_'.join(col) for col in run_s.columns]

    rid = get_rid(f)
    run_s["rid"] = rid
    return run_s.reset_index()

def eval_proportion_to_pop(df_stats, pop = "p1"):
    df = df_stats.sort_values(by = ["level_0", "focal_pop_id", "selection_bins"])
    df.index = range(df.shape[0])

    d0 = df[df.focal_pop_id == pop]
    d0.index = range(d0.shape[0])
    dd0 = pd.concat([d0] * df.focal_pop_id.unique().shape[0])
    dd0 = dd0.sort_values(by = ["level_0", "focal_pop_id", "selection_bins"])
    dd0.index = range(dd0.shape[0])


    ratios = df.iloc[:,3:] / dd0.iloc[:,3:]
    ratios = pd.concat([df.iloc[:,:3],
                        ratios], 
                       axis = 1)

    return ratios

def get_rid(f):
    return int(f.split("-")[-1].split('.')[0])



if __name__ == "__main__":


    path = '../data/'
    list_files = os.listdir(path)

    fns = [f for f in list_files if 'l-1000000_' in f]
    fh = [f for f in fns if 'huber' in f]
    f = fh[3]

    %time summaries = pd.concat([get_simul_summaries(f) for f in fns])

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

    x = 
    m_dt.dominance.isinf()
x.hist()

x = pd.DataFrame() 
x["selection"] = [-1.0, -0.1, -0.01, -0.001, -0.0001, -0.00000]
x

m_dt[m_dt.dominance < 0].loc[:,["selection", "dominance"]]

huber_dominance(x)


