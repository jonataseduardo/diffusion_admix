#Effects of Selection on Recent and Strong Admixture of Human Populations

In this document we will analyse how strong and recent admixture events effects
the allele frequency spectrum, genetic load and other related summary statists. 

Here we will show the results of simulations of humans Out-Of-Africa expansion
the European population inferred of by Jouganous2017. In this model, humans first
expand its effective population size from 11273 individuals to 23721 in the
African Continent around 302 kyears (1 generation is equal 29 years). The
Out-of-Africa split takes place 125 kyears in the past with a population with
3004 effective size. The European population starts to exponentially increase
in its size around 42 kyears after starting from initial size of 2271 when a
bottleneck took place. Finally, we create 3 admixed population from a single
pulse admixture between the African and European population for with 0.25, 0.5,
0.75 proportions of Africa ancestry.  

## WF simulations with recombination and DFE (SLIM3)

In theses simulations we used used the a gamma distributed DFE with shape =
-0.0065 and rate = 0.186. The mutations rate was set to 1e-17 and the
recombination rate used was 1e-8.  

In the following figures the African poulation label is p1 the European
population label is p2, p3 and p4 have 25\%, 50\%, 75\%  of Africa ancestry.  

The boxplots sumarizes the results of 100 simulations. 


### Morton Genetic load  

The Morton genetic load is defined as the average number of deleterious
mutations by individual in a population wich can be computed from as the sum of
the mean frequency of deleterious mutation (or sum of the first frequency
moment). 

[]!(mu_1_sum_bs.png)

### Genetic load 

The genetic load is defined as L = -s sum \[2h x1 + (1 - 2h) x2] 


### Genetic diversity

The genetic diversity is defined as the summation of the nucleotide
heterosigosity 


### Number of segregating sites 
