import numpy as np

def simple_gamma(mu = 1e-6, 
                 rho = 1e-6, 
                 dominance = 0.5,
                 alpha = -0.00657402, 
                 beta = 0.186,
                 lenght = 1000,
                 **kwargs
                 ):

    init = """ 
        initialize() {{
            initializemutationrate({mu});
            initializerecombinationrate({rho});
            initializemutationtype("m1", {dominance}, "g", {alpha}, {beta});
            initializegenomicelementtype("g1", m1, 1.0);
            initializegenomicelement(g1, 0, {lenght});
        }}
        """
    init = init.format(mu = mu, rho = rho, dominance = dominance,
                       alpha = alpha, beta = beta, lenght = lenght)

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
        
    if dominance_modifier is "huber_model": 
        init = init + huber_model(**kwargs)

    return init
    


if __name__ == "__main__":
    init = init_mutations() 
    print init
