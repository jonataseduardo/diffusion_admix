#!/bin/python3

import numpy

numpy.set_printoptions(suppress=True)

NSAMPLE = [100]
GAMMA = - numpy.round(numpy.logspace(-2, 1.5, 15), decimals = 3)
DOMINANCE = numpy.round(numpy.linspace(0, 0.5, 11), decimals = 3)

LENGHT = 1000000

afr_eur_moments = expand('data/AfEu_ns{nsample}_g{gamma}_d{dominance}.fs', 
                          nsample=NSAMPLE,
                          gamma=GAMMA,
                          dominance=DOMINANCE)

afr_eur_admix_huber = expand('data/AfEuAdmx_huber_l-{lenght}_rid-{rid}.txt',
                              lenght = LENGHT, 
                              rid = range(100))

afr_eur_admix = expand('data/AfEuAdmx_l-{lenght}_rid-{rid}.txt',
                        lenght = LENGHT, 
                        rid = range(100))


afr_eur_admix = expand('data/AfEuAdmx_d-{dominance}_l-{lenght}_rid-{rid}.txt',
                        lenght = LENGHT, 
                        dominance = [0.0, 0.05, 0.2],
                        rid = range(100))

rule all:
    input:
        afr_eur_moments, 
        afr_eur_admix_huber,
        afr_eur_admix

rule run_afr_eur_moments:
    output:
        'data/AfEu_ns{nsample}_g{gamma}_d{dominance}.fs'
    shell:
        "python2.7 -W ignore"
        " scripts/AfrEur_split.py" 
        " -ns {wildcards.nsample}" 
        " -g {wildcards.gamma}"
        " -d {wildcards.dominance}"
        " -o {output}"

rule run_afr_eur_admix_huber:
    output:
        'data/AfEuAdmx_huber_l-{lenght}_rid-{rid}.txt'
    shell:
        "python"
        " scripts/slim_simuls.py" 
        " -g jouganous2017"
        " -mm simple_gamma"
        " -dm huber_model" 
        " -l {wildcards.lenght}"
        " -o {output}"
        
rule run_afr_eur_admix:
    output:
        'data/AfEuAdmx_l-{lenght}_rid-{rid}.txt'
    shell:
        "python"
        " scripts/slim_simuls.py" 
        " -g jouganous2017"
        " -mm simple_gamma"
        " -l {wildcards.lenght}"
        " -o {output}"

rule run_afr_eur_admix_dominance:
    output:
        'data/AfEuAdmx_d-{dominance}_l-{lenght}_rid-{rid}.txt'
    shell:
        "python"
        " scripts/slim_simuls.py" 
        " -g jouganous2017"
        " -mm simple_gamma"
        " -d {wildcards.dominance}"
        " -l {wildcards.lenght}"
        " -o {output}"
