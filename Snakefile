#!/bin/python3

import numpy

numpy.set_printoptions(suppress=True)

<<<<<<< HEAD
NSAMPLE = [100]
#GAMMA = - numpy.round(numpy.logspace(-2, 1.5, 15), decimals = 3)
GAMMA = - numpy.round(numpy.logspace(-1, 3, 30), decimals = 3)
DOMINANCE_M = numpy.round(numpy.linspace(0, 0.5, 11), decimals = 3)
DOMINANCE_S = numpy.round(numpy.linspace(0, 0.5, 11), decimals = 3)
LENGHT = 1000000

afr_eur_moments = expand('moments_data/AfEu_ns{nsample}_g{gamma}_d{dominance}.fs', 
                          nsample=NSAMPLE,
                          gamma=GAMMA,
                          dominance=DOMINANCE_M)

afr_eur_admix_huber = expand('slim_data/AfEuAdmx_huber_l-{lenght}_rid-{rid}.txt',
                              lenght = LENGHT, 
                              rid = range(100))

afr_eur_admix_n = expand('slim_data/AfEuAdmx_neutral_l-{lenght}_rid-{rid}.txt',
                          lenght = LENGHT, 
                          rid = range(100))

afr_eur_admix_d = expand('slim_data/AfEuAdmx_d-{dominance}_l-{lenght}_rid-{rid}.txt',
                          lenght = LENGHT, 
                          dominance = DOMINANCE_S,
                          rid = range(100))

rule all:
    input:
        #afr_eur_moments, 
        #afr_eur_admix_huber,
        afr_eur_admix_n,
        afr_eur_admix_d

=======
>>>>>>> e8f72771c9a041f1ef3910d481bcf53e47f844c1

rule run_afr_eur_moments:
    output:
        'moments_data/AfEu_ns{nsample}_g{gamma}_d{dominance}.fs'
    shell:
        "python2.7 -W ignore"
        " scripts/AfrEur_split.py" 
        " -ns {wildcards.nsample}" 
        " -g {wildcards.gamma}"
        " -d {wildcards.dominance}"
        " -o {output}"

rule run_afr_eur_admix_huber:
    output:
        'slim_data/AfEuAdmx_huber_l-{lenght}_rid-{rid}.txt'
    shell:
        "python"
        " scripts/slim_simuls.py" 
        " -g jouganous2017"
        " -mm simple_gamma"
        " -dm huber_model" 
        " -l {wildcards.lenght}"
        " -o {output}"
        
rule run_afr_eur_admix_neutral:
    output:
        'slim_data/AfEuAdmx_neutral_l-{lenght}_rid-{rid}.txt'
    shell:
        "python"
        " scripts/slim_simuls.py" 
        " -g jouganous2017"
        " -mm neutral"
        " -l {wildcards.lenght}"
        " -o {output}"

rule run_afr_eur_admix_dominance:
    output:
        'slim_data/AfEuAdmx_d-{dominance}_l-{lenght}_rid-{rid}.txt'
    shell:
        "python"
        " scripts/slim_simuls.py" 
        " -g jouganous2017"
        " -mm simple_gamma"
        " -d {wildcards.dominance}"
        " -l {wildcards.lenght}"
        " -o {output}"
