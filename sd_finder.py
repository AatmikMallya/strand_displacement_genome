#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from tqdm.notebook import tqdm
# from ipywidgets import IntProgress
from IPython.display import display
import time
import datetime
import itertools
from statistics import mean, stdev
import rapidfuzz
import string
import os.path
import string
from lru import LRU
tab = str.maketrans("ACTG", "TGAC")


# In[55]:


l = LRU(1024+256+64+16+4) 


# In[53]:


l['A']='T'
l['C']='G'
l['T']='A'
l['G']='C'
l['AA']='TT'


# In[54]:


l['TT']='AA'


# In[56]:


sequences = ['A','T','C','G']
for length in range(1,6):
    for seq in itertools.product(sequences, repeat=length):
        seq_str = ''.join(seq)
        l[seq_str] = seq_str.translate(tab)[::-1]


# In[58]:


len(l)


# In[ ]:





# In[ ]:





# In[ ]:





# ## Helper functions for SD finder

# In[2]:


# Read sequences from a fasta file into an array
def load_fasta(file):
    fasta_seqs = SeqIO.parse(open(file), 'fasta')
    seqs = [seq.seq for seq in fasta_seqs]
    return seqs


# In[3]:


# finds first index at which substring occurs within string, allowing for mismatches
def find_fuzzy(string, substring, threshold):
    str_len = len(substring)    
    for i in range(len(string) - str_len + 1):
        sub_str = string[i:i+str_len]
        if rapidfuzz.distance.Hamming.distance(sub_str, substring) <= threshold: 
            return i
    return -1


# In[4]:


# Same as find_fuzzy but also adds a normalized hamming distance constraint on the substring
def find_fuzzy_extra(string, substring, threshold, p_distance):
    str_len = len(substring)    
    for i in range(len(string) - str_len + 1):
        sub_str = string[i:i+str_len]
        if rapidfuzz.distance.Hamming.distance(sub_str, substring) <= threshold and rapidfuzz.distance.Hamming.normalized_distance(sub_str, substring) <= p_distance: 
            return i
    return -1


# In[5]:


def a_astar_a(dna):
    longest_substring = ""
    longest_length = 0
    
    for i in range(len(dna)):
        for j in range(i, len(dna)):
            # Extract the substring
            substring = dna[i:j+1]
            reverse_complement = substring.reverse_complement()
            # Check if the reverse complement substring exists after the substring
            index_rev_comp = dna.find(reverse_complement, j+1)

            # Check if the substring occurs again after the reverse complement
            if index_rev_comp != -1 and dna.find(substring, index_rev_comp+len(reverse_complement)) != -1:
                # Update the longest substring and its length if necessary
                if j-i+1 > longest_length:
                    longest_substring = substring
                    longest_length = j-i+1
    return len(longest_substring)


# In[6]:


def astar_a_a(dna):
    longest_substring = ""
    longest_length = 0
    
    for i in range(len(dna)):
        for j in range(i, len(dna)):
            # Extract the substring
            substring = dna[i:j+1]
            reverse_complement = substring.reverse_complement()
            # Check if the reverse complement substring exists after the substring
            index_rev_comp = dna.find(reverse_complement, j+1)

            # Check if the substring occurs again after the reverse complement
            if index_rev_comp != -1 and dna.find(reverse_complement, index_rev_comp+len(reverse_complement)) != -1:
                # Update the longest substring and its length if necessary
                if j-i+1 > longest_length:
                    longest_substring = substring
                    longest_length = j-i+1
    return len(longest_substring)


# In[7]:


def a_a_astar(dna):
    longest_substring = ""
    longest_length = 0
    
    for i in range(len(dna)):
        for j in range(i, len(dna)):
            # Extract the substring
            substring = dna[i:j+1]
            
            # Check if the reverse complement substring exists after the substring
            index_2nd_substr = dna.find(substring, j+1)
            
            reverse_complement = substring.reverse_complement()

            # Check if the substring occurs again after the reverse complement
            if index_2nd_substr != -1 and dna.find(reverse_complement, index_2nd_substr+len(reverse_complement)) != -1:
                # Update the longest substring and its length if necessary
                if j-i+1 > longest_length:
                    longest_substring = substring
                    longest_length = j-i+1
    return len(longest_substring)


# In[59]:


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
    return len(longest_substring)


# In[60]:


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
    return len(longest_substring)


# In[61]:


def a_a_astar_threshold(dna, threshold, p_distance):
    length = astar_a_a_threshold2(dna[::-1],threshold,p_distance)
    return length
    # return length, seq.complement()


# In[76]:


def a_astar_a_threshold2(dna, threshold, p_distance):
    longest_substring = ""
    longest_length = 0
    
    for i in range(len(dna)):
        for j in range(i, len(dna)):
            # Extract the substring
            substring = dna[i:j+1]
            try: 
                reverse_complement = l[substring]
            except KeyError:
                reverse_complement = substring.translate(tab)[::-1]
            # Check if the reverse complement substring exists after the substring
            index_rev_comp = dna.find(reverse_complement, j+1)

            # Check if the substring occurs again after the reverse complement
            if index_rev_comp != -1 and find_fuzzy_extra(dna[index_rev_comp+len(reverse_complement):], substring, threshold, p_distance) != -1:
                # Update the longest substring and its length if necessary
                if j-i+1 > longest_length:
                    longest_substring = substring
                    longest_length = j-i+1
    return len(longest_substring)


# In[ ]:





# In[ ]:





# ## Benchmarking

# In[181]:


timeit a_astar_a(Seq(''.join(np.random.choice(['A','C','T','G'], 200))))


# In[182]:


timeit a_a_astar(Seq(''.join(np.random.choice(['A','C','T','G'], 200))))


# In[183]:


timeit astar_a_a(Seq(''.join(np.random.choice(['A','C','T','G'], 200))))


# In[83]:


timeit a_astar_a_threshold(''.join(np.random.choice(['A','C','T','G'], 300)), 3, 0.5)


# In[72]:


timeit a_a_astar_threshold(''.join(np.random.choice(['A','C','T','G'], 200)), 3, 0.5)


# In[73]:


timeit astar_a_a_threshold(''.join(np.random.choice(['A','C','T','G'], 200)), 3, 0.5)


# In[ ]:





# In[79]:


timeit a_astar_a_threshold2(''.join(np.random.choice(['A','C','T','G'], 200)), 3, 0.5)


# In[27]:


timeit a_a_astar_threshold2(''.join(np.random.choice(['A','C','T','G'], 200)), 3, 0.5)


# In[28]:


timeit astar_a_a_threshold2(''.join(np.random.choice(['A','C','T','G'], 200)), 3, 0.5)


# ## Data sources

# In[20]:


# RF00162: SAM riboswitch https://rfam.org/family/RF00162
# RF03072: raiA riboswitch https://rfam.org/family/RF03072
# RF01739: glnA riboswitch https://rfam.org/family/RF01739
# RF03054: Xanthine riboswitch https://rfam.org/family/RF03054
# RF02913: pemK cis-reg https://rfam.org/family/RF02913
# RF01750: ZMP/ZTP riboswitch
# RF01852: tRNA-Sec https://rfam.org/family/RF01852
# RF01854: Bacteria_large_SRP https://rfam.org/family/RF01854
# RF02541: LSU_rRNA_bacteria https://rfam.org/family/RF02541
# RF00001: 5S_rRNA https://rfam.org/family/RF00001
# RF00005: tRNA https://rfam.org/family/RF00005

rfam_codes = ['RF00162','RF03072','RF01739','RF03054','RF02913','RF01750','RF01852','RF01854','RF02541','RF00001','RF00005']


# ## Main functions

# In[12]:


# For each sequence in the fasta file, finds the longest stretch out of the three categories and saves it in a df
def run_fasta(file, outfile=False):
    start = time.time()
    seqs = load_fasta(file)

    lengths = []
    seq_lengths = []
    for seq in seqs:
        seq = seq.replace('U','T')
        if not set(seq).issubset({'A','C','T','G'}) or len(seq) > 300:
            continue

        max_len = max(a_astar_a(seq), astar_a_a(seq), a_a_astar(seq))
        
        lengths.append(max_len)
        seq_lengths.append(len(seq))

    data = {'length':lengths, 'seq_length':seq_lengths}
    df = pd.DataFrame.from_dict(data)
    
    if outfile:
        df.to_csv(outfile, index=False)
    
    print(f'Total time: {time.time() - start} seconds.')
    return df


# In[13]:


# For each sequence in the fasta file, finds the longest stretch out of the three categories and saves it in a df
# threshold = number of mismatches that are allowed between the invading strand and its complimentary strand (hamming distance)
# p_distance = percentage of mismatches that are allowed btwn invading and complimentary strand (norm hamming distance)
def run_fasta_threshold(file, threshold, p_distance, outfile=False):
    start = time.time()
    seqs = load_fasta(file)

    lengths = []
    seq_lengths = []
    for seq in seqs:
        seq = seq.replace('U','T')
        if not set(seq).issubset({'A','C','T','G'}) or len(seq) > 300:
            continue
            
        max1 = a_astar_a_threshold(seq, threshold, p_distance)
        max2 = astar_a_a_threshold(seq, threshold, p_distance)
        max3 = a_a_astar_threshold(seq, threshold, p_distance)
        max_len = max(max1, max2, max3)
        
        lengths.append(max_len)
        seq_lengths.append(len(seq))

    data = {'length':lengths, 'seq_length':seq_lengths}
    df = pd.DataFrame.from_dict(data)
    
    if outfile:
        df.to_csv(outfile, index=False)
    
    print(f'Total time: {time.time() - start} seconds.')
    # return df


# In[14]:


# Run script on random sequences
def run_random_threshold(threshold, p_distance, outfile):
    seq_len_range = [i for i in range(25, 351, 25)]

    lengths = []
    seq_lengths = []
    begin = time.time()
    for seq_len in seq_len_range:
        start = time.time()
        num_trials = 100
        for _ in range(num_trials):
            seq = Seq(''.join(random.choices('ACTG',k=seq_len)))
            
            max1 = a_astar_a_threshold(seq, threshold, p_distance)
            max2 = astar_a_a_threshold(seq, threshold, p_distance)
            max3 = a_a_astar_threshold(seq, threshold, p_distance)
            max_len = max(max1, max2, max3)
            
            lengths.append(max_len)
            seq_lengths.append(seq_len)

        end = time.time()

        # print(f'Length {seq_len} complete. Time: {datetime.timedelta(seconds=end-start)}. Total time elapsed: {datetime.timedelta(seconds=end-begin)}')

    data = {'length':lengths, 'seq_length':seq_lengths}
    df = pd.DataFrame.from_dict(data)

    df.to_csv(outfile, index=False)


# In[ ]:





# In[22]:


for threshold in range(6):
    print('On threshold: ' + str(threshold))
    for code in rfam_codes:
        rfam_file = f'results/{code}_threshold_{threshold}_extra.csv'
        if not os.path.exists(rfam_file):
            run_fasta_threshold(f'{code}.fa', threshold, 0.5, outfile=rfam_file)
    
    random_seq_file = f'results/random_threshold_{threshold}_extra.csv'
    if not os.path.exists(random_seq_file):
        run_random_threshold(threshold, 0.5, random_seq_file)


# In[ ]:





# In[ ]:


run_fasta_threshold('RF03072.fa', 2, outfile=f'RF03072_threshold_0.csv')


# In[135]:


run_fasta('RF03072.fa', outfile='RF03072_test.csv')


# In[ ]:





# In[8]:


def plot_comparison(filename, label):
    df_rfam = pd.read_csv(filename)
    df_random = pd.read_csv('random_seqs.csv')
    # For smoothing the plot
    df_rfam['seq_length'] = df_rfam['seq_length'].apply(lambda x: round(x, -1))

    plt.title('SD detection')
    plt.xlabel('Length of sequence')
    plt.ylabel('Length of longest stretch')
    sns.lineplot(data=df_rfam, x='seq_length', y='length', label=label)
    sns.lineplot(data=df_random, x='seq_length', y='length', label='random sequences')


# In[9]:


plot_comparison('RF00162.csv', 'SAM riboswitch (RF00162)')


# In[10]:


plot_comparison('RF01739.csv', 'glnA riboswitch (RF01739)')


# In[11]:


plot_comparison('RF02913.csv', 'raiA riboswitch (RF02913)')


# In[12]:


plot_comparison('RF01750.csv', 'ZMP/ZTP riboswitch (RF01750)')


# In[13]:


plot_comparison('RF01739.csv', 'pemK cis-reg (RF01739)')


# In[193]:


plot_comparison('RF03072.csv', 'tRNA-Sec (RF03072)')


# In[14]:


# Run script on random sequences
def run_random():
    seq_len_range = [i for i in range(25, 351, 25)]

    lengths = []
    seq_lengths = []
    begin = time.time()
    for seq_len in seq_len_range:
        start = time.time()
        num_trials = 100
        for _ in range(num_trials):
            seq = Seq(''.join(random.choices('ACTG',k=seq_len)))
            max_len = max(a_astar_a(seq), astar_a_a(seq), a_a_astar(seq))
            lengths.append(max_len)
            seq_lengths.append(seq_len)

        end = time.time()

        print(f'Length {seq_len} complete. Time: {datetime.timedelta(seconds=end-start)}. Total time elapsed: {datetime.timedelta(seconds=end-begin)}')

    data = {'length':lengths, 'seq_length':seq_lengths}
    df = pd.DataFrame.from_dict(data)

df.to_csv('random_seqs.csv')


# In[166]:





# In[15]:





# In[153]:


df = pd.read_csv('RF00162_threshold_5.csv')


# In[154]:


df.sort_values(by='seq_length')


# In[141]:


seqs = load_fasta('RF00162.fa')
seqs = [seq.replace('U','T') for seq in seqs]
seqs = [seq for seq in seqs if set(seq).issubset({'A','C','T','G'})]


# In[94]:


seq = seqs[4815]


# In[140]:


str(seq)


# In[137]:


a_star_i = find_fuzzy(seq[seq.find(a.reverse_complement())+len(a.reverse_complement()):], a, 10) + seq.find(a.reverse_complement())+len(a.reverse_complement())
a_star = seq[a_star_i:a_star_i+len(a)]


# In[130]:


a, length = test(seq, 10)
str(a), length


# In[138]:


print(a + '\n' + a_star)


# In[ ]:





# In[117]:


seq.find(a), seq.find(a.reverse_complement())


# In[113]:


matrix = build_bp_matrix(seq, False, False)


# In[114]:


plt.figure(figsize=(19.2, 14.4))
title = 'Base pair matrix'
plt.title(title)
plt.imshow(matrix, cmap='hot')
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




