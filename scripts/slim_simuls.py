from subprocess import Popen, PIPE
import tempfile, os

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

def balick_split_model(N0 = 1000, N1 = 100, Tburn = 100, 
                       Tbn = 100, Tpost_bn = 100): 

    demography = """
        1 {{ sim.addSubpop("p1", {N0}); }}
        {T1}{{sim.addSubpopSplit("p2", {N1}, p1); }}
        {T2}{{p2.setSubpopulationSize({N0}); }}
        {T3}{{
            p1.outputSample({N0}); 
            p2.outputSample({N0}); 
            }}
        """
    return demography.format(N0 = N0, N1 = N1, 
                             T1 = Tburn, 
                             T2 = Tburn + Tbn, 
                             T3 = Tburn + Tbn + Tpost_bn)

def run_slim(slim_cmd, data_path = None):
    out = None
    with tempfile.NamedTemporaryFile() as file:
        file.write(slim_cmd)
        file.delete = False
    try: 
        process =  Popen(["slim", file.name], 
                         stdout = PIPE, stderr = PIPE, stdin = PIPE)
        out = process.communicate()
    finally:
        os.remove(file.name)

    return out

if __name__ == "__main__":
    init = init_mutations(mutation_model = "neutral", length = 1000, dominance = 0.3)
    demography = balick_split_model()
    print init + demography


    x = run_slim(init + demography)
    print x[0]
    x.communicate()
