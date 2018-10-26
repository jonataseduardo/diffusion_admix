import numpy
import pandas
import moments
from itertools import combinations

def OutOfAfrica_stepwise(s, 
                         h,
                         (n1, n2, n3),
                         (stepsAf, stepsB, stepsEuAs)):
    """
    A three-population model used to model out-of-Africa demography.
    """

    params = np.array([2.10065897, 0.25066579, 0.22247642, 3.05297944,
                       0.09022469, 5.82773903, 3.79104318, 0.25730946,
                       0.12569788, 1.07182332, 0.36429414, 0.1108222, 
                       0.07072507])

    (nuAf, nuB, nuEu0, nuEu, nuAs0, nuAs, 
     mAfB, mAfEu, mAfAs, mEuAs, TAf, TB, TEuAs) = params

    Ttotal = TAf + TB + TEuAs

    dtAf = TAf / stepsAf
    dtB = TB / stepsB
    dtEuAs = TEuAs / stepsEuAs

    fs_history = OrderedDict()

    #first step: a single population
    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2+n3)
    fs = moments.Spectrum(sts)
    
    dt = dtAf 
    time = 0
    step = 0 
    split_1 = split_2 = True
    while(time <= 1.001 * Ttotal ):
        time += dt

        if(time <= 1.001 * TAf):
            print('before split' , step, time )
            #integrate for time TAf (with constant population)
            fs.integrate([nuAf], dt, min(0.5 * dt, 0.05), gamma = s[0] )
            fs_history[step] = [time, fs]
            step += 1

        if(time > TAf and time <= 1.001 * (TAf + TB)):
            if(split_1):
                print('first split' , step, time )
                split_1 = False
                dt = dtB
                fs = moments.Manips.split_1D_to_2D(fs, n1, n2+n3)
                fs_history[step] = [time, fs]
                step += 1
                mig1=numpy.array([[0, mAfB],[mAfB, 0]])

            print('after first split' , step, time )
            fs.integrate([nuAf, nuB], dt, min(0.5 * dt, 0.05), m=mig1, 
                    gamma = s[0:2])
            fs_history[step] = [time, fs]
            step += 1

        if(time > TAf + TB and time <= 1.001 * (TAf + TB + TEuAs)):
            if(split_2):
                print('sencond split' , step, time )
                split_2 = False
                dt = dtEuAs
                fs = moments.Manips.split_2D_to_3D_2(fs, n2, n3)
                fs_history[step] = [time, fs]
                step += 1
                # migration rates matrix
                mig2=numpy.array([[0, mAfEu, mAfAs],
                                  [mAfEu, 0, mEuAs],
                                  [mAfAs, mEuAs, 0]])
                nuEu0t, nuAs0t = nuEu0, nuAs0

            print('after sencond split' , step, time )

            #define functions for population sizes
            nuEu_func = lambda t: nuEu0t*(nuEu/nuEu0t)**(t/TEuAs)
            nuAs_func = lambda t: nuAs0t*(nuAs/nuAs0t)**(t/TEuAs)
            nu2 = lambda t: [nuAf, nuEu_func(t), nuAs_func(t)]

            # migration rates matrix
            mig2=numpy.array([[0, mAfEu, mAfAs],[mAfEu, 0, mEuAs],[mAfAs, mEuAs, 0]])
            fs.integrate(nu2, dt, min(0.5 * dt, 0.05), m=mig2, 
                         gamma = s)

            fs_history[step] = [time, fs]
            nuEu0t = nuEu0t*(nuEu/nuEu0t)**(dt/TEuAs)
            nuAs0t = nuAs0t*(nuAs/nuAs0t)**(dt/TEuAs)
            step += 1
                                
    return fs_history

def OOA_stats(fs_history, s, h):

    stats_pop = dict()
    stats_pop["step"] = []
    stats_pop["time"] = []
    stats_pop["pop_id"] = []
    stats_pop["pi"] = []
    stats_pop["load"] = []
    stats_pop["efficacy"] = []

    for step in fs_history.keys():
        np = fs_history[step][1].Npop
        rnp = range(np)
        lpop = [(list(set(rnp) ^ set(i)), i) for i in combinations(rnp, np - 1)]
        for i in lpop:
            stats_pop["step"] += [step]
            stats_pop["time"] += [fs_history[step][0]]
            stats_pop["pop_id"] += [i[0][0]]
            try:
                fs_i = fs_history[step][1].marginalize(i[1])
            except:
                fs_i = fs_history[step][1]

            stats_pop["pi"] += [fs_i.pi()]
            stats_pop["load"] += [mutation_load(fs_i, s[i[0][0]], h[i[0][0]])]
            stats_pop["efficacy"] += [efficacy_of_selection(fs_i, s[i[0][0]], h[i[0][0]])]

    return pandas.DataFrame.from_dict(stats_pop)
