#!/usr/bin/env python
# coding: utf-8


import numpy as np
import random
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns

# Read sequences from a fasta file into an array
def load_fasta(file):
    fasta_seqs = SeqIO.parse(open(file), 'fasta')
    seqs = [str(seq.seq) for seq in fasta_seqs]
    return seqs

# Check if two base pairs are complementary (wobble include G-T interactions)
def are_compatible(a, b, include_wobble=False):
    if a == 'U':
        a = 'T'
    if b == 'U':
        b = 'T'
    
    if ((a == 'C' and b == 'G') or (a == 'G' and b == 'C') or 
        (a == 'A' and b == 'T') or (a == 'T' and b == 'A')):
            return True

    if include_wobble and ((a == 'G' and b == 'T') or (a == 'T' and b == 'G')):
        return True

    return False

# Finds the longest stretch of two reverse-complementary sub-sequences in a sequence
def find_longest_stretch(seq):
    n = len(seq)
    matrix = np.zeros((n,n))
    maxlen = 0
    res_i = 0
    res_j = 0
    include_wobble = False

    # Populate matrix with base pairs
    for i in range(n):
        for j in range(i,n):
            if are_compatible(seq[i], seq[j], include_wobble):
                matrix[i,j] = 1
                
    # Find longest stretch in matrix
    for i in range(n):
        for j in range(i,n):
            sublen = 0
            k = 0
            while i-k>=0 and k+j<n:
                if matrix[i-k,j+k] == 1:
                    sublen += 1
                    k += 1
                else:
                    break
            if sublen > maxlen:
                maxlen = sublen
                res_i = i
                res_j = j

    a_start = res_i - maxlen + 1
    a_end = res_i
    b_start = res_j
    b_end = res_j + maxlen - 1

    indices = [a_start, a_end, b_start, b_end]
    return matrix, indices

# Show a heatmap of the base pairs, with the longest stretch highlighted
def show_matrix(matrix, indices):
    a_start, a_end, b_start, b_end = indices
    print(f'Sequence A: {a_start},{a_end}')
    print(f'Sequence B: {b_start},{b_end}')
    
    for i in range(a_end - a_start + 1):
        matrix[a_start+i, b_end-i] = 0.4

    plt.figure(figsize=(19.2, 14.4))
    plt.title('Base pair matrix (red = longest reverse complementary stretch)')
    plt.imshow(matrix, cmap='hot')
    plt.show()

def main():
    seqs = load_fasta('random_seq.fasta')
    matrix, indices = find_longest_stretch(seqs[0])
    show_matrix(matrix, indices)
main()

