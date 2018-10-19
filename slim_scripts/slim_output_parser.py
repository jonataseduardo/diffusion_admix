"""
File: slim_output_parser.py
Author: Jonatas Cesar   
Email: jonataseduardo@gmail.com
Github: https://github.com/jonataseduardo   
Description: Parse the default output SliM3.
For a detailed description of the output see 
section 4.2 of the SliM manual
"""
import pandas
import numpy

def slim_output_parser(file_name):
    """
    Description: Parse the default output SliM3.
    For a detailed description of the output see 
    section 4.2 of the SliM manual

    param file_name: path to the Slim3 output file
    ptype: string
    return (mutations_df, genomes_df, individuals_df):
    rtype: (pandas.DataFrame, pandas.DataFrame, pandas.DataFrame)
    """

    with open(file_name) as f:
        content = f.readlines()

    content = [line.rstrip() for line in content]

    line_key = None
    slim_data = dict()
    slim_data["Mutations:"] = []
    slim_data["Individuals:"] = []
    slim_data["Genomes:"] = []

    pop_id = 0
    for line in content:
        if line.startswith("#OUT"):
            lsp =  line.split()
            pop_id = lsp[3]
            sample_size = lsp[4]
        if line == "Mutations:":
           line_key = line 
        if line == "Individuals:":
           line_key = line 
        if line == "Genomes:":
           line_key = line 
        if line_key is not None:
            if line_key != line and not line.startswith("#"): 
                slim_data[line_key] += [line.split() + [pop_id, sample_size]]
     
    mutations_df = pandas.DataFrame(slim_data["Mutations:"])
    genomes_df = pandas.DataFrame(slim_data["Genomes:"])
    individuals_df = pandas.DataFrame(slim_data["Individuals:"])

    m_columns = {0: 'prm_mid',
                 1: 'tmp_mid',       
                 2: 'm_type',        
                 3: 'position',      
                 4: 'selection',     
                 5: 'dominance',     
                 6: 'origin_pop',    
                 7: 'birth_time',    
                 8: 'allele_count',
                 9: 'focal_pop_id',
                 10: 'sample_size'}

    mutations_df.rename(columns = m_columns, inplace = True )

    def change_to_numeric(x):
        try:
            z =  pandas.to_numeric(x)
        except:
            z = x
        return z

    mutations_df = mutations_df.apply(change_to_numeric)

    return (mutations_df, genomes_df, individuals_df)

def make_slim_sfs(mutations_df):
    """
    Evaluate the Site Freqrency Spectrum of Slim3 simulations. 

    Params mutation_df: pandas.DataFrame with position focal_pop_id and
    allele_count colums

    Return numpy.array: site frequency spectrum of populations 
    """

    #Eventually the more the one mutation can appear in one loci. 
    #see slim documentation introduction
    #aggfunc = max selects the oldest mutation with high probability
    aux = mutations_df.pivot_table(index = ['position'],  
                                   columns = 'focal_pop_id',
                                   values = 'allele_count', 
                                   aggfunc = max ).fillna(0)

    cols = ['p' + str(i + 1) for i in range(aux.shape[1])]

    sfs_vals = aux.groupby(cols).agg('size')
    sfs_vals.name = 'num_counts'

    sfs_vals = sfs_vals.reset_index()

    #the pivot_table is converting an integer to a flot.
    sfs_vals = sfs_vals.apply(lambda x : pandas.to_numeric(x, downcast = 'signed'))

    dim = mutations_df.groupby('focal_pop_id').agg({'sample_size': max})

    sfs = numpy.zeros(dim.values.reshape(-1) + 1, dtype = int)

    idx = [sfs_vals.loc[:, i].values for i in cols]

    sfs[tuple(idx)] = sfs_vals.num_counts.values

    return (sfs, sfs_vals)
