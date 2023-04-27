#!/usr/bin/env python
# coding: utf-8

import numpy as np
import random
from BCBio.GFF import GFFExaminer
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
# import matplotlib.pyplot as plt
# import seaborn as sns
import pandas as pd
from tqdm.notebook import tqdm
import time
import datetime
import itertools
from statistics import mean, stdev
import rapidfuzz
import string
import os.path
import string
tab = str.maketrans("ACTG", "TGAC")


## Helper functions for SD finder
# Read sequences from a fasta file into an array
def load_fasta(file):
    fasta_seqs = SeqIO.parse(open(file), 'fasta')
    seqs = [str(seq.seq) for seq in fasta_seqs]
    return seqs

# Same as find_fuzzy but also adds a normalized hamming distance constraint on the substring
def find_fuzzy_extra(string, substring, threshold, p_distance):
    str_len = len(substring)    
    for i in range(len(string) - str_len + 1):
        sub_str = string[i:i+str_len]
        if rapidfuzz.distance.Hamming.distance(sub_str, substring) <= threshold and rapidfuzz.distance.Hamming.normalized_distance(sub_str, substring) <= p_distance: 
            return i
    return -1

def a_astar_a_threshold(dna, threshold, p_distance):
    longest_substring = ""
    longest_length = 0
    
    for i in range(len(dna)):
        for j in range(i, len(dna)):
            # Extract the substring
            substring = dna[i:j+1]
            reverse_complement = substring.translate(tab)[::-1]
            # Check if the reverse complement substring exists after the substring
            index_rev_comp = dna.find(reverse_complement, j+1)

            # Check if the substring occurs again after the reverse complement
            if index_rev_comp != -1 and find_fuzzy_extra(dna[index_rev_comp+len(reverse_complement):], substring, threshold, p_distance) != -1:
                # Update the longest substring and its length if necessary
                if j-i+1 > longest_length:
                    longest_substring = substring
                    longest_length = j-i+1
    return len(longest_substring), longest_substring

def astar_a_a_threshold(dna, threshold, p_distance):
    longest_substring = ""
    longest_length = 0
    
    for i in range(len(dna)):
        for j in range(i, len(dna)):
            # Extract the substring
            substring = dna[i:j+1]
            reverse_complement = substring.translate(tab)[::-1]
            # Check if the reverse complement substring exists after the substring
            index_rev_comp = dna.find(reverse_complement, j+1)

            # Check if the substring occurs again after the reverse complement
            
            if index_rev_comp != -1 and find_fuzzy_extra(dna[index_rev_comp+len(reverse_complement):], reverse_complement, threshold, p_distance) != -1:
                # Update the longest substring and its length if necessary
                if j-i+1 > longest_length:
                    longest_substring = substring
                    longest_length = j-i+1
    return len(longest_substring), longest_substring

def a_a_astar_threshold(dna, threshold, p_distance):
    length, seq = astar_a_a_threshold(dna[::-1],threshold,p_distance)
    return length, seq.translate(tab)


# Only keep standard bases (A,C,T,G) in the genome
def preprocess_genome(genome, window_size):
    indices = [i for i, base in enumerate(genome) if base not in ['A','C','T','G']]
    if len(indices) == 1:
        i = indices[0]
        new_genome = genome[:i-window_size] + genome[i-1+window_size:]
        return new_genome
    else:
        raise Exception('TODO')

# Run strand displacement script over a genome
def run_SD_genome(genome, stride, window_size, threshold, p_distance, outfile):
    window_size = 150
    # Split genome into equal segments of length 'window_size'
    # seqs = [genome[i:i+window_size] for i in range(0,len(genome),window_size)]

    lengths = []
    sd_seqs = []
    window_idx = []
    
    start = time.time()
    for i in range(len(genome) - window_size + 1):
        seq = genome[i:i+window_size]
        if i % 100 == 0:
            print(i)

        max1, seq1 = a_astar_a_threshold(seq, threshold, p_distance)
        max2, seq2 = astar_a_a_threshold(seq, threshold, p_distance)
        max3, seq3 = a_a_astar_threshold(seq, threshold, p_distance)

        max_len = max(max1, max2, max3)
        if max_len == max1:
            sd_seq = seq1
        elif max_len == max2:
            sd_seq = seq2
        else:
            sd_seq = seq3
        
        lengths.append(max_len)
        sd_seqs.append(str(sd_seq))
       

    data = {'length':lengths, 'sd_seq':sd_seqs}
    df = pd.DataFrame.from_dict(data)

    df.to_csv(outfile, index=False)
    
    end = time.time()
    print(f'Total run time: {end - start} seconds')

# Splits the genome into equal sections for parallel processing on cluster (optional)
def select_genome_split(genome, num_splits, split_idx):
    # Determine the length of each split
    split_len = len(genome) // num_splits
    
    # Determine the starting and ending indices of the selected split
    start_idx = (split_idx - 1) * split_len
    end_idx = start_idx + split_len
    
    # If this is the last split, include any remaining characters
    if split_idx == num_splits:
        end_idx = len(genome)
    
    selected_split = genome[start_idx:end_idx]
    
    return selected_split


# ## Bacteria

# http://bacteria.ensembl.org/Escherichia_coli_w_gca_000184185/Info/Index
ecoli_genome = '/home/amallya2/strand_displacement/data/Escherichia_coli_w_gca_000184185.ASM18418v1.dna.chromosome.Chromosome.fa'
genome = load_fasta(ecoli_genome)[0]

genome = load_fasta(ecoli_genome)[0]
stride = 150
window_size = 150
threshold = 5
p_distance = 0.5
outfile = '/home/amallya2/strand_displacement/ecoli_sd.csv'


genome = load_fasta(ecoli_genome)[0]
indices = [i for i, base in enumerate(genome) if base not in ['A','C','T','G']]


ecoli_gff = "Escherichia_coli_w_gca_000184185.ASM18418v1.56.chromosome.Chromosome.gff3"
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
df_ribo = df2[df2.attributes.str.contains('riboswitch')]
df_ribo['window_size'] = df_ribo.end - df_ribo.start



lengths = []
sd_seqs = []
for i, row in df_ribo.iterrows():
    seq = genome[row.start:row.end]

    max1, seq1 = a_astar_a_threshold(seq, threshold, p_distance)
    max2, seq2 = astar_a_a_threshold(seq, threshold, p_distance)
    max3, seq3 = a_a_astar_threshold(seq, threshold, p_distance)

    max_len = max(max1, max2, max3)
    if max_len == max1:
        sd_seq = seq1
    elif max_len == max2:
        sd_seq = seq2
    else:
        sd_seq = seq3

    lengths.append(max_len)
    sd_seqs.append(str(sd_seq))

df_result = pd.DataFrame.from_dict({'length':lengths, 'sd_seq':sd_seqs, })

