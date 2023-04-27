import numpy as np
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pandasql as ps
import time
import datetime
import itertools
from statistics import mean, stdev
import rapidfuzz
import string
import pprint
from BCBio.GFF import GFFExaminer
from BCBio import GFF


dfs = []
for i in range(1,9):
    df = pd.read_csv(f'ecoli_250window_{i}.csv')
    df['start_idx'] = df.index + (612371+250) * (i-1)
    dfs.append(df)
    
    
df = pd.concat(dfs)


ecoli_gff = "data/Escherichia_coli_w_gca_000184185.ASM18418v1.56.chromosome.Chromosome.gff3"
# parse the GFF file using BCBio
with open(ecoli_gff) as handle:
    records = list(GFF.parse(handle))

# create an empty dataframe to store the results
df_list = []

# iterate over each record in the parsed GFF file
for record in records:
    for feature in record.features:
        # get the start and end coordinates of the feature
        start = feature.location.start.position
        end = feature.location.end.position

        # get the strand of the feature
        if feature.location.strand == 1:
            strand = '+'
        elif feature.location.strand == -1:
            strand = '-'
        else:
            strand = '.'

        # get the feature type and attributes
        feature_type = feature.type
        attributes = feature.qualifiers

        # get the sequence ID from the record
        seqid = record.id

        # create a new dataframe row with the extracted data
        row = pd.DataFrame({'seqid': [seqid], 'start': [start], 'end': [end], 
                            'strand': [strand], 'feature': [feature_type], 
                            'attributes': [attributes]})

        # append the row dataframe to the list
        df_list.append(row)

# concatenate all the dataframes in the list into a single dataframe
df2 = pd.concat(df_list, ignore_index=True)
df2['attributes'] = df2['attributes'].apply(str)
df2 = df2.drop(11)



query = '''
    SELECT *
    FROM df
    INNER JOIN df2
    ON df.start_idx BETWEEN df2.start AND df2.end
'''

# execute query with pandasql
result_df = ps.sqldf(query)

result_df.to_csv('combine_gff_results2.csv')











