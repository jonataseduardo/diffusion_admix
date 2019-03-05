#!/bin/python3

import numpy

numpy.set_printoptions(suppress=True)

NSAMPLE = [5]
#GAMMA = [-100, -50, -10, -5, -1, -0.5, -0.1, -0.05, -0.01]
#DOMINANCE = [0.0, 0.01, 0.05, 0.1, 0.3, 0.5]
GAMMA = - numpy.round(numpy.logspace(-2, 2, 2), decimals = 2)
DOMINANCE = numpy.round(numpy.linspace(0, 0.5, 2), decimals = 2)

#rule all:
#    input:
#        expand('data/AfEu_ns{nsample}_g{gamma}_d{dominance}.fs', 
#               nsample=NSAMPLE,
#               gamma=GAMMA,
#               dominance=DOMINANCE)

rule moments_afr_eur:
    input:
        expand('data/AfEu_ns{nsample}_g{gamma}_d{dominance}.fs', 
               nsample=NSAMPLE,
               gamma=GAMMA,
               dominance=DOMINANCE)
    output:
        'data/AfEu_ns{nsample}_g{gamma}_d{dominance}.fs'
    shell:
        "python2.7 -W ignore"
        " scripts/AfrEur_split.py" 
        " -ns {wildcards.nsample}" 
        " -g {wildcards.gamma}"
        " -d {wildcards.dominance}"
        " -o {output}"
