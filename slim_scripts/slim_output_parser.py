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

    def change_to_numeric(x):
        try:
            z =  pandas.to_numeric(x)
        except:
            z = x
        return z

    mutations_df = mutations_df.apply(change_to_numeric)

    return (mutations_df, genomes_df, individuals_df)

mutations_df = slim_output_parser(file_name)[0] 
mutations_df = slim_output_parser("split2.txt")[0] 

mutations_df[mutations_df.position == 10]


def make_slim_sfs(mutations_df):

    mutations_df.groupby(["focal_pop_id", "allele_count"]
                         ).agg({"position":"count"}).reset_index()


    aux = mutations_df.pivot_table(index = ['position'],  
                                   columns = 'focal_pop_id',
                                   values = 'allele_count', 
                                   aggfunc = sum ).fillna(0)
    aux.shape
    aux.head()

    cols = ['p' + str(i + 1) for i in range(aux.shape[1])]
    cols

    sfs_vals = aux.groupby(cols).agg('size')

    sfs_vals.reset_index()

    mutations_df.groupby(["allele_count"]).agg({"position":"count"})



