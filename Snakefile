#!/bin/python3

import numpy

numpy.set_printoptions(suppress=True)


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
        
rule run_afr_eur_admix_neutral:
    output:
        'data/AfEuAdmx_neutral_l-{lenght}_rid-{rid}.txt'
    shell:
        "python"
        " scripts/slim_simuls.py" 
        " -g jouganous2017"
        " -mm neutral"
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
