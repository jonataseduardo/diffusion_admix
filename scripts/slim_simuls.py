from subprocess import Popen, PIPE
import tempfile, os
import click

def neutral(mu = 1e-7, 
            rho = 1e-8, 
            length = 1000,
            **kwargs
            ):

    init = """ 
        initialize() {{
            initializeMutationRate({mu});
            initializeRecombinationRate({rho});
            initializeMutationType("m1", 0.5, "f", 0.0);
            initializeGenomicElementType("g1", m1, 1.0);
            initializeGenomicElement(g1, 0, {length});
        }}
        """
    init = init.format(mu = mu, rho = rho, length = length)

    return init

def simple_gamma(mu = 1e-7, 
                 rho = 1e-8, 
                 dominance = 0.5,
                 alpha = -0.00657402, 
                 beta = 0.186,
                 length = 1000,
                 **kwargs
                 ):

    init = """ 
        initialize() {{
            initializeMutationRate({mu});
            initializeRecombinationRate({rho});
            initializeMutationType("m1", {dominance}, "g", {alpha}, {beta});
            initializeGenomicElementType("g1", m1, 1.0);
            initializeGenomicElement(g1, 0, {length});
        }}
        """
    init = init.format(mu = mu, rho = rho, dominance = dominance,
                       alpha = alpha, beta = beta, length = length)

    return init

def synthetic_chr(mu = 1e-7, 
                  rho = 1e-8, 
                  alpha = -0.00657402, 
                  beta = 0.186,
                  length = 1000000,
                  **kwargs
                  ):
     init = """ 
        initialize() {{
            initializeMutationRate({mu});

            // deleterious
            initializeMutationType("m1", 0.5, "g", {alpha}, {beta}); 
            initializeMutationType("m2", 0.5, "f", 0.0);         // synonymous
            initializeMutationType("m3", 0.5, "f", 0.0);         // non-coding
            initializeMutationType("m4", 0.5, "e", 0.1);         // beneficial
            
            initializeGenomicElementType("g1", c(m2,m1,m4), c(2,8,0.1)); // exon
            initializeGenomicElementType("g2", c(m3,m1), c(9,1));      // intron
            initializeGenomicElementType("g3", c(m3), 1);          // non-coding
            
            // Generate random genes along an approx {length}-base chr
            base = 0;
            
            while (base < {length}) {{
                // make a non-coding region
                nc_length = asInteger(runif(1, 100, 5000));
                initializeGenomicElement(g3, base, base + nc_length - 1);
                base = base + nc_length;
                
                // make first exon
                ex_length = asInteger(rlnorm(1, log(50), log(2))) + 1;
                initializeGenomicElement(g1, base, base + ex_length - 1);
                base = base + ex_length;
                
                // make additional intron-exon pairs
                do
                {{
                    in_length = asInteger(rlnorm(1, log(100), log(1.5))) + 10;
                    initializeGenomicElement(g2, base, base + in_length - 1);
                    base = base + in_length;
                    
                    ex_length = asInteger(rlnorm(1, log(50), log(2))) + 1;
                    initializeGenomicElement(g1, base, base + ex_length - 1);
                    base = base + ex_length;
                }}
                while (runif(1) < 0.8);  // 20% probability of stopping
            }}
            
            // final non-coding region
            nc_length = asInteger(runif(1, 100, 5000));
            initializeGenomicElement(g3, base, base + nc_length - 1);
            
            // single recombination rate
            initializeRecombinationRate({rho});
        }}
        """
     init = init.format(mu = mu, rho = rho, alpha = alpha, beta = beta, 
                        length = length)
     return init

def huber_model(h_intercept = 0.5, h_rate = 1e6, **kwargs):

    hs_modifier = """
        fitness(m1) {{
          if (homozygous)
            return 1.0 + mut.selectionCoeff;
          else
            return 1.0 + mut.selectionCoeff /
                   ( {intercept} + {h_rate} * mut.selectionCoeff);
        }}
        """
    return hs_modifier.format(intercept = 1. / h_intercept, h_rate = h_rate)

def out_of_africa(n1 = 216, n2 = 198, n3 = 206, **kwargs):
    demographics = """
        // Create the ancestral African population
        1 {{ sim.addSubpop("p1", 7310); }}

        // Expand the African population to 14474
        // This occurs 148000 years (5920) generations ago
        52080 {{ p1.setSubpopulationSize(14474); }}

        // Split non-Africans from Africans and set up migration between them
        // This occurs 51000 years (2040 generations) ago
        55960 {{
                sim.addSubpopSplit("p2", 1861, p1);
                p1.setMigrationRates(c(p2), c(15e-5));
                p2.setMigrationRates(c(p1), c(15e-5));
        }}

        // Split p2 into European and East Asian subpopulations
        // This occurs 23000 years (920 generations) ago
        57080 {{
                sim.addSubpopSplit("p3", 554, p2);
                p2.setSubpopulationSize(1032);  // reduce European size

                // Set migration rates for the rest of the simulation
                p1.setMigrationRates(c(p2, p3), c(2.5e-5, 0.78e-5));
                p2.setMigrationRates(c(p1, p3), c(2.5e-5, 3.11e-5));
                p3.setMigrationRates(c(p1, p2), c(0.78e-5, 3.11e-5));
        }}

        // Set up exponential growth in Europe and East Asia
        // Where N(0) is the base subpopulation size and t = gen - 57080:
        //    N(Europe) should be int(round(N(0) * e^(0.0038*t)))
        //    N(East Asia) should be int(round(N(0) * e^(0.0048*t)))
        57080:58000 {{
                t = sim.generation - 57080;
                p2_size = round(1032 * exp(0.0038 * t));
                p3_size = round(554 * exp(0.0048 * t));
                
                p2.setSubpopulationSize(asInteger(p2_size));
                p3.setSubpopulationSize(asInteger(p3_size));
        }}

        // Generation 58000 is the present.  Output and terminate.
        58000 late() {{
                p1.outputSample({n1}); // YRI phase 3 sample of size 108
                p2.outputSample({n2}); // CEU phase 3 sample of size 99
                p3.outputSample({n3}); // CHB phase 3 sample of size 103
        }}
        """
    return demographics.format(n1 = n1, n2 = n2, n3 = n3)

def afr_eur_admix(n = 1000, **kwargs):
    demographics = """
        // Create the ancestral African population
        1 {{ sim.addSubpop("p1", 11273); }}

        // Expand the African population to 23721
        // This occurs 312000 years (10758) generations ago
        47242 {{ p1.setSubpopulationSize(23721); }}

        // Split non-Africans from Africans and set up migration between them
        // This occurs 125 kyears (4310 generations) ago
        53690 {{
                sim.addSubpopSplit("p2", 3104, p1);
                p1.setMigrationRates(c(p2), c(15.8e-5));
                p2.setMigrationRates(c(p1), c(15.8e-5));
        }}

        // European subpopulation second bottleneck 
        // This occurs 42.3 kyears (1459 generations) ago
        56541 {{
                p2.setSubpopulationSize(2271);  // reduce European size

                // Set migration rates for the rest of the simulation
                p1.setMigrationRates(c(p2), c(1.1e-5));
                p2.setMigrationRates(c(p1), c(1.1e-5));
        }}

        // Set up exponential growth in Europe and East Asia
        // Where N(0) is the base subpopulation size and t = gen - 56541:
        //    N(Europe) should be int(round(N(0) * e^(0.00196*t)))
        56541:58000 {{
                t = sim.generation - 56541;
                p2_size = round(2271 * exp(0.00196 * t));
                p2.setSubpopulationSize(asInteger(p2_size));
        }}

        58001{{
                sim.addSubpop("p3", {n});
                p3.setMigrationRates(c(p1,p2), c(0.25, 0.75));
                sim.addSubpop("p4", {n});
                p4.setMigrationRates(c(p1,p2), c(0.5, 0.5));
                sim.addSubpop("p5", {n});
                p5.setMigrationRates(c(p1,p2), c(0.75, 0.25));
              }}

        // Generation 58001 is the present. Output and terminate.
        58001 late() {{
                p1.outputSample({n}); 
                p2.outputSample({n}); 
                p3.outputSample({n}); 
                p4.outputSample({n}); 
                p5.outputSample({n}); 
        }}
        """
    return demographics.format(n = n)

def balick_split_model(N0 = 1000, N1 = 100, 
                       Tburn = 100, Tbn = 100, Tpost_bn = 100, 
                       **kwargs): 

    demography = """
        1 {{ sim.addSubpop("p1", {N0}); }}
        {T1}{{sim.addSubpopSplit("p2", {N1}, p1); }}
        {T2}{{p2.setSubpopulationSize({N0}); }}
        {T3}{{sim.addSubpop("p3", {N0});
              p3.setMigrationRates(c(p1,p2), c(0.25, 0.75));
              sim.addSubpop("p4", {N0});
              p4.setMigrationRates(c(p1,p2), c(0.5, 0.5));
              sim.addSubpop("p5", {N0});
              p5.setMigrationRates(c(p1,p2), c(0.75, 0.25));
              }}

        {T3} late() {{
            p1.outputSample({N0}); 
            p2.outputSample({N0}); 
            p3.outputSample({N0}); 
            p4.outputSample({N0}); 
            p5.outputSample({N0}); 
            }}
        """
    return demography.format(N0 = N0, N1 = N1, 
                             T1 = Tburn, 
                             T2 = Tburn + Tbn, 
                             T3 = Tburn + Tbn + Tpost_bn)

def single_pop(N0 = 1000, Tburn = 100):

    demography = """
        1 {{ sim.addSubpop("p1", {N0}); }}

        {T1} late() {{
            p1.outputSample({N0}); 
            }}
        """
    return demography.format(N0 = N0, T1 = Tburn)

def run_slim(slim_cmd, output = None):
    out = None
    with tempfile.NamedTemporaryFile(delete = False) as file:
        file.write(slim_cmd.encode())
    try: 
        if output is None:
            process =  Popen(["slim", file.name], 
                             stdout = PIPE, stderr = PIPE)
        else:
            cmd = "slim " + file.name + " > " + output
            process =  Popen(cmd, 
                             stdout = PIPE, stderr = PIPE, shell = True)
        out = process.communicate()
    finally:
        os.remove(file.name)

    return out

@click.command()
@click.option("-mm", "--mutation-model", 
              help = "Mutation Model",
              default="neutral", 
              show_default=True,
              required=True, 
              type=click.Choice(["neutral", "simple_gamma", "synthetic_chr"])
              )

@click.option("-mu", "--mutation-rate", 
              help = "Mutation rate per pair base",
              default=1.44e-8, 
             show_default=True,
              required=True, 
              type=float
              )

@click.option("-rho", "--recombination-rate", 
              help = "recombination rate per pair base",
              default=1e-8, 
              show_default=True,
              required=True, 
              type=float
              )

@click.option("-l", "--length", 
              help = "number of pair bases to be simulated",
              default=int(1e4), 
              show_default=True,
              required=True, 
              type=int
              )

@click.option("-a", "--alpha", 
              help = "rate of the gamma distribution",
              default=-0.00657402, 
              show_default=True,
              required=True, 
              type=float
              )

@click.option("-b", "--beta", 
              help = "shape of the gamma distribution",
              default=0.186, 
              show_default=True,
              required=True, 
              type=float
              )

@click.option("-dm", "--dominance-model", 
              help = "Dominance model",
              default="constant", 
              show_default=True,
              required=False, 
              type=click.Choice(["constant", "huber_model", "logistic"])
              )

@click.option("-d", "--dominance-value", 
              help = "Dominance value considering constant dominance model",
              default= 0.5, 
              show_default=True,
              required=False, 
              type=float
              )

@click.option("-g", "--demographic-model", 
              help = "Demographic Models",
              default="single_pop", 
              show_default=True,
              required=False, 
              type=click.Choice(["single_pop", 
                                 "balick2015", 
                                 "gravel2013", 
                                 "jouganous2017"]
                                )
              )
@click.option("-N0", "--pop0-size", 
              help = "Ancestral pop size for single_pop or balick2015 models",
              default=int(1e3), 
              show_default=True,
              required=True, 
              type=int
              )
@click.option("-N1", "--pop1-size", 
              help = "bottleneck pop size balick2015 model",
              default=int(1e2), 
              show_default=True,
              required=True, 
              type=int
              )
@click.option("-t0", "--time-burn", 
              help = "burn time for single_pop or balick2015 models",
              default=int(1e3), 
              show_default=True,
              required=True, 
              type=int
              )
@click.option("-tb", "--time-bottleneck", 
              help = "bottleneck time in the balick2015 model",
              default=int(1e3), 
              show_default=True,
              required=True, 
              type=int
              )
@click.option("-tp", "--time-postbottleneck", 
              help = "post bottleneck time the balick2015 model", 
              default=int(1e3), 
              show_default=True,
              required=True, 
              type=int
              )
@click.option("-o", "--output", 
              help = "output file name", 
              default= "output.txt",
              show_default=True,
              type=str
              )
def main(mutation_model,
         dominance_model,
         mutation_rate, 
         recombination_rate,
         length,
         dominance_value,
         alpha,
         beta,
         demographic_model,
         pop0_size,
         pop1_size,
         time_burn,
         time_bottleneck,
         time_postbottleneck,
         output,
         **kwargs
         ):

    if mutation_model == "simple_gamma":
        init = simple_gamma(mu = mutation_rate, 
                            rho = recombination_rate, 
                            length = length, 
                            dominance = dominance_value, 
                            alpha = alpha, 
                            beta = beta,
                            **kwargs)

    elif mutation_model == "neutral":
        init = neutral(mu = mutation_rate, 
                       rho = recombination_rate, 
                       length = length, 
                       **kwargs)

    elif mutation_model == "synthetic_chr":
        init = synthetic_chr(mu = mutation_rate, 
                             rho = recombination_rate, 
                             length = length, 
                             dominance = dominance_value,
                             alpha = alpha, 
                             beta = beta,
                             **kwargs)

    if dominance_model == "huber_model": 
        init = init + huber_model(**kwargs)

    if demographic_model == "gravel2013":
        demography = out_of_africa()
    elif demographic_model == "jouganous2017":
        demography = afr_eur_admix()
    elif demographic_model == "single_pop":
        demography = single_pop(N0 = pop0_size,
                                Tburn = time_burn)
    elif demographic_model == "balick2017":
        demography = balick_split_model(N0 = pop0_size,
                                        N1 = pop1_size, 
                                        Tburn = time_burn,
                                        Tbn = time_bottleneck,
                                        Tpost_bn = time_postbottleneck
                                        )
                                        
    simul_code = init + demography
    #click.echo(simul_code)
                                 
    slim_error = run_slim(simul_code, output = output)[1]
    if slim_error != '':
        click.echo(slim_error)
    pass


if __name__ == "__main__":

    main()
