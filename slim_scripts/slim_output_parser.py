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
        if line == "Mutations:":
           line_key = line 
           pop_id += 1
        if line == "Individuals:":
           line_key = line 
        if line == "Genomes:":
           line_key = line 
        if line_key is not None:
            if line_key != line and not line.startswith("#"): 
                slim_data[line_key] += [line.split() + ['p' + str(pop_id)]]
     
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
                 9: 'focal_pop_id'}

    mutations_df.rename(columns = m_columns, inplace = True )

    return (mutations_df, genomes_df, individuals_df)
