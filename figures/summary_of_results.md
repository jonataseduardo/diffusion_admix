#Results of SLIM3 simulations

In this document, we will look to some effects of negative selection under
different dominance constraints on the genetic variability of human admixed
populations. For this, we used the humans Out-Of-Africa demography inferred of
by Jouganous2017 adding single pulse admixture events in the end. 

In this model, humans first expand their effective population size from 11273
individuals to 23721 in the African Continent around 302 kyears (1 generation
is equal 29 years). The Out-of-Africa split takes place 125 kyears in the past
with a population with 3004 effective size. The European population starts to
grow in size exponentially around 42kyears. At the end of the demographic
processes, we create 3 admixed population from a single pulse admixture between
the African and European population with 0.25, 0.5, 0.75 African proportion.

Under this demography, we simulated chromosomes of 1Mbp with deleterious
mutations with selection coefficient sampled from a gamma distribution DFE
(shape = -0.0065 and rate = 0.186). The mutations rate was set to 1e-7 and the
recombination rate used was 1e-8. We have run simulations with several values
of dominance, varying from h = 0 (totally recessive) to h = 0.5 (additive)
wherein each simulation the value of the dominance coefficient is constant. We
also simulated the model proposed by Huber et.al. 2018 where the dominance is a
function of selection h = 1/(a + b \* s) where h(s=0) = 1/a = 0.5 is the value
of dominance of neutral mutations and b = 1e6. 

In the following figures, the African population label is p1 the European
population label is p2. The labels p3, p4, and p5 correspond populations
with 25\%, 50\%, 75\%  of African ancestry.  

The boxplots summarize the results of 100 simulations, and in each simulation,
we sampled 1000 individuals from each population. The quadrants in the figures
shows the resulting statistics evaluated on sets of deleterious variants
aggregated by the magnitude of the selection coefficient. 


### Morton Genetic load  

The Morton genetic load is defined as the average number of deleterious
mutations by individual in a population which can be computed from as the sum
of the mean frequency of deleterious mutation (or sum of the first frequency
moment). 

![](slim-mu1-bs.png)

### Genetic load 

The genetic load is defined as L = -s sum \[2h x1 + (1 - 2h) x2] 
(I forgot to multiply the load value by -1 in the figure)

![](load_sum_bs.png)

### Genetic diversity

The genetic diversity is defined as the summation of the nucleotide
heterozygousity 

![](htz_sum_bs.png)

### Number of segregating sites 

![](selection_count_bs.png)


### Morton efficacy

The Morton measures the rate in which selection removes deleterious mutations
in a populations. 

![](morton_sum_bs.png)

## Conclusions 

* It is possible to see an accumulation of extremely deleterious mutations in
  the simulated European population in comparison with the African population.
  This accumulation is very pronounced when the mutations are totally
  recessive, but it is still present for values of h up to 0.2 (I needs to
  check this statement more carefully).  

* Most of the statistics for the admixed populations are between the values for
  the European and African populations. 

* The number of segregating sites in the admixed populations is higher then in
  the parental populations. 
  

## Caveats 

* The DFE used in this simulations is narrower then the one estimated for
  humans in the recent literature. 
* The coefficient b used also is too big making most variants almost recessive.
* The simulation takes too much time to run, probably I will need to rescale
  the effective pop size so it can run faster. 

##  Next steps

* Replicate this results on the gnomAD data set and check if there is an
  accumulation of deleterious mutation for European and East Asian population
  using the rank of score of several functional annotations such as CADD, Gerp,
  etc as a proxy for selection coefficients. 

* Replicate the results of the slim simulations using WF diffusion simulations.
  This will need to be done using dadi since the Moments get numeric errors for
  values of 2Ns < -50 

#Gnomad Data

## Phenotypic Scores 

![](./normalized_scores_quantililes.png)

## Mean Frequency of deleterious mutations

![](./mean_AF_afr_nfe.png)

## SFS projection

![](./mean_project_nfe.png)

![](./numsites_project_nfe.png)

![](./var_project_nfe.png)





