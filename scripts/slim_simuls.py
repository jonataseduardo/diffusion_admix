
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
                  dominance = 0.5, 
                  alpha = -0.00657402, 
                  beta = 0.186,
                  length = 1000000,
                  **kwargs
                  ):

    init = """ 
        initialize() {{
            initializeMutationRate({mu});

            // deleterious
            initializeMutationType("m1", {dominance}, "g", {alpha}, {beta}); 
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
    init = init.format(mu = mu, rho = rho, dominance = dominance,
                       alpha = alpha, beta = beta, length = length)

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

def init_mutations(mutation_model = "simple_gamma", 
                   dominance_modifier = None,
                   **kwargs
                   ):

    if mutation_model is "simple_gamma":
        init = simple_gamma(**kwargs)
    elif mutation_model is "neutral":
        init = neutral(**kwargs)
    elif mutation_model is "synthetic_chr":
        init = synthetic_chr(**kwargs)
        
    if dominance_modifier is "huber_model": 
        init = init + huber_model(**kwargs)

    return init

def out_of_africa(

        ):
        """
        // Create the ancestral African population
        1 { sim.addSubpop("p1", 7310); }

        // Expand the African population to 14474
        // This occurs 148000 years (5920) generations ago
        52080 { p1.setSubpopulationSize(14474); }

        // Split non-Africans from Africans and set up migration between them
        // This occurs 51000 years (2040 generations) ago
        55960 {
                sim.addSubpopSplit("p2", 1861, p1);
                p1.setMigrationRates(c(p2), c(15e-5));
                p2.setMigrationRates(c(p1), c(15e-5));
        }

        // Split p2 into European and East Asian subpopulations
        // This occurs 23000 years (920 generations) ago
        57080 {
                sim.addSubpopSplit("p3", 554, p2);
                p2.setSubpopulationSize(1032);  // reduce European size

                // Set migration rates for the rest of the simulation
                p1.setMigrationRates(c(p2, p3), c(2.5e-5, 0.78e-5));
                p2.setMigrationRates(c(p1, p3), c(2.5e-5, 3.11e-5));
                p3.setMigrationRates(c(p1, p2), c(0.78e-5, 3.11e-5));
        }

        // Set up exponential growth in Europe and East Asia
        // Where N(0) is the base subpopulation size and t = gen - 57080:
        //    N(Europe) should be int(round(N(0) * e^(0.0038*t)))
        //    N(East Asia) should be int(round(N(0) * e^(0.0048*t)))
        57080:58000 {
                t = sim.generation - 57080;
                p2_size = round(1032 * exp(0.0038 * t));
                p3_size = round(554 * exp(0.0048 * t));
                
                p2.setSubpopulationSize(asInteger(p2_size));
                p3.setSubpopulationSize(asInteger(p3_size));
        }

        // Generation 58000 is the present.  Output and terminate.
        //58000 late() {
        //	p1.outputSample(216); // YRI phase 3 sample of size 108
        //	p2.outputSample(198); // CEU phase 3 sample of size 99
        //	p3.outputSample(206); // CHB phase 3 sample of size 103
        //}

        // Generation 58000 is the present.  Output and terminate.
        58000 late() {
                p1.outputSample(6); // YRI phase 3 sample of size 108
                p2.outputSample(8); // CEU phase 3 sample of size 99
                p3.outputSample(6); // CHB phase 3 sample of size 103
        }
        """

if __name__ == "__main__":
    init = init_mutations(mutation_model = "synthetic_chr", dominance = 0.3)
    print init
