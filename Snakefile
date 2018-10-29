#!/bin/python3

NSAMPLE = [100]
GAMMA = [-100, -50, -10, -5, -1, -0.5, -0.1, -0.05, -0.01]
DOMINANCE = [0.0, 0.01, 0.05, 0.1, 0.3, 0.5]
FADMIX = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

rule all:
    input:
        expand('data/AFW_ns{nsample}_g{gamma}_d{dominance}_f{fadmix}.fs', 
               nsample=NSAMPLE,
               gamma=GAMMA,
               dominance=DOMINANCE,
               fadmix=FADMIX)

rule moments_simul:
    output:
        'data/AFW_ns{nsample}_g{gamma}_d{dominance}_f{fadmix}.fs'
    shell:
        "python2.7 -W ignore"
        " scripts/AfrEur_admix.py" 
        " -ns {wildcards.nsample}" 
        " -g {wildcards.gamma}"
        " -d {wildcards.dominance}"
        " -f {wildcards.fadmix}"
        " -o {output}"
